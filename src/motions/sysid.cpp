#include "motions/sysid.hpp"
#include "chassis.hpp"
#include "main.h"
#include <cmath>
#include <cstdio>
#include <random>

SystemIdentification::SystemIdentification(Chassis &chassis)
    : chassis_(chassis), fp(nullptr), _currentTime(0.0) {}

void SystemIdentification::runSession(bool isAngular, const char *filename) {
  if (!beginLog(filename))
    return;

  double uL_scale = 12.0;
  double uR_scale = isAngular ? -12.0 : 12.0;

  // 1. Static Phase (Measure sensor noise floor)
  printf("  [Phase 1/4]  Static logging (Don't move!)...\n");
  for (int i = 0; i < (int)(2.0 / DT); ++i) {
    auto row = sampleRow(_currentTime, 0, 0);
    logRow(row, "static");
    _currentTime += DT;
    pros::delay((int)(DT * 1000));
  }

  // 2. MultiStep (DC Gains)
  printf("  [Phase 2/4]  Running Multi-step...\n");
  runMultiStep(uL_scale, uR_scale, 0.1, 10);

  // Stop segment for cleaner fitting segments
  for (int i = 0; i < 10; ++i) {
    auto row = sampleRow(_currentTime, 0, 0);
    logRow(row, "stop");
    _currentTime += DT;
    pros::delay((int)(DT * 1000));
  }
  std::cout << "finished" << std::endl;
  waitForButton();

  // 3. PRBS (High-freq motor dynamics, softened)
  printf("  [Phase 3/4]  Running PRBS...\n");
  runPRBS(uL_scale * 0.7, uR_scale * 0.7, 0.15, 10.0);

  // Stop segment
  for (int i = 0; i < 10; ++i) {
    auto row = sampleRow(_currentTime, 0, 0);
    logRow(row, "stop");
    _currentTime += DT;
    pros::delay((int)(DT * 1000));
  }
  waitForButton();

  // 4. Chirp (Medium-freq dynamics)
  printf("  [Phase 4/4]  Running Chirp...\n");
  runChirp(uL_scale * 0.8, uR_scale * 0.8, 0.1, 5.0, 5.0);

  endLog();
  printf("SysID session complete: %s\n", filename);
}

bool SystemIdentification::beginLog(const char *filename) {
  fp = fopen(filename, "w");
  if (!fp) {
    printf("ERROR: Could not open %s\n", filename);
    return false;
  }
  fprintf(fp, "t,uL,uR,uL_act,uR_act,vL,vR,heading,tag\n");
  return true;
}

void SystemIdentification::endLog() {
  if (fp)
    fclose(fp);
  fp = nullptr;
  chassis_.tank(0, 0);
}

void SystemIdentification::logRow(const RawRow &row, const char *tag) {
  if (!fp)
    return;
  fprintf(fp, "%.4f,%.2f,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%s\n", row.t, row.uL,
          row.uR, row.u_act_l, row.u_act_r, row.vL, row.vR, row.heading, tag);
}

SystemIdentification::RawRow
SystemIdentification::sampleRow(double t, double uL, double uR) {
  RawRow row;
  row.t = t;
  row.uL = uL;
  row.uR = uR;

  // Actual voltage averages
  double uL_sum = 0, uR_sum = 0;
  int cl = 0, cr = 0;
  for (auto &m : chassis_.left_motors) {
    double v = m->get_voltage();
    if (v != PROS_ERR_F) {
      uL_sum += v;
      cl++;
    }
  }
  for (auto &m : chassis_.right_motors) {
    double v = m->get_voltage();
    if (v != PROS_ERR_F) {
      uR_sum += v;
      cr++;
    }
  }
  row.u_act_l = (cl > 0 ? (uL_sum / cl) / 1000.0 : 0.0);
  row.u_act_r = (cr > 0 ? (uR_sum / cr) / 1000.0 : 0.0);

  row.vL = chassis_.get_raw_left_velocity();
  row.vR = chassis_.get_raw_right_velocity();
  row.heading = chassis_.get_rotation();
  return row;
}

void SystemIdentification::runMultiStep(double maxUL, double maxUR,
                                        double stepSize, int numSteps) {
  for (int step = 0; step < numSteps; ++step) {
    double uL = (step + 1) * stepSize * maxUL;
    double uR = (step + 1) * stepSize * maxUR;
    for (int i = 0; i < (int)(1.0 / DT); ++i) {
      chassis_.tank((int)(uL * 1000.0), (int)(uR * 1000.0));
      auto row = sampleRow(_currentTime, uL, uR);
      logRow(row, "multistep");
      _currentTime += DT;
      pros::delay((int)(DT * 1000));
    }
    chassis_.tank(0, 0);
    for (int i = 0; i < (int)(0.5 / DT); ++i) {
      auto row = sampleRow(_currentTime, 0, 0);
      logRow(row, "coastdown");
      _currentTime += DT;
      pros::delay((int)(DT * 1000));
    }
    for (int i = 0; i < 5; ++i) {
      auto row = sampleRow(_currentTime, 0, 0);
      logRow(row, "stop");
      _currentTime += DT;
      pros::delay((int)(DT * 1000));
    }
    waitForButton();
  }
}

void SystemIdentification::runPRBS(double maxUL, double maxUR, double holdTime,
                                   double duration) {
  double uL = 0, uR = 0;
  uint32_t lastSwitchTime = 0;
  int mode = 0;

  for (int i = 0; i < (int)(duration / DT); ++i) {
    if (pros::millis() - lastSwitchTime > holdTime * 1000) {
      // Cycle through modes to ensure zero-mean movement (stay in place)
      switch (mode % 4) {
      case 0:
        uL = maxUL;
        uR = maxUR;
        break; // Forward
      case 1:
        uL = -maxUL;
        uR = -maxUR;
        break; // Backward
      case 2:
        uL = maxUL;
        uR = -maxUR;
        break; // Turn Right
      case 3:
        uL = -maxUL;
        uR = maxUR;
        break; // Turn Left
      }
      mode++;
      lastSwitchTime = pros::millis();
    }
    chassis_.tank((int)(uL * 1000.0), (int)(uR * 1000.0));
    auto row = sampleRow(_currentTime, uL, uR);
    logRow(row, "prbs");
    _currentTime += DT;
    pros::delay((int)(DT * 1000));
  }
  chassis_.tank(0, 0);
}

void SystemIdentification::runChirp(double maxUL, double maxUR, double minFreq,
                                    double maxFreq, double duration) {
  for (int i = 0; i < (int)(duration / DT); ++i) {
    double local_t = i * DT;
    double freq = minFreq + (maxFreq - minFreq) * (local_t / duration);
    double uL = maxUL * std::sin(2.0 * M_PI * freq * local_t);
    double uR = maxUR * std::sin(2.0 * M_PI * freq * local_t);
    chassis_.tank((int)(uL * 1000.0), (int)(uR * 1000.0));
    auto row = sampleRow(_currentTime, uL, uR);
    logRow(row, "chirp");
    _currentTime += DT;
    pros::delay((int)(DT * 1000));
  }
}

void SystemIdentification::waitForButton() {
  pros::Controller master(pros::E_CONTROLLER_MASTER);
  printf("  [Step Complete]  Press (A) to continue...\n");
  while (!master.get_digital(pros::E_CONTROLLER_DIGITAL_A)) {
    pros::delay(20);
  }
  while (master.get_digital(pros::E_CONTROLLER_DIGITAL_A)) {
    pros::delay(20);
  }
  pros::delay(400); // Set-up time
}

void run_linear_sysid(Chassis &chassis) {
  SystemIdentification si(chassis);
  si.runSession(false, "/usd/sysid_linear.csv");
  chassis.tank(0, 0);
}

void run_angular_sysid(Chassis &chassis) {
  SystemIdentification si(chassis);
  si.runSession(true, "/usd/sysid_angular.csv");
  chassis.tank(0, 0);
}
