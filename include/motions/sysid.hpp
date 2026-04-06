#ifndef ACCURATE_GAIN_MEASURE
#define ACCURATE_GAIN_MEASURE

#include "chassis.hpp"
#include <cstdio>
#include <vector>

/**
 * @brief System Identification — data collection and logging.
 *
 * Runs PRBS + chirp + multi-step excitation on both linear and angular axes,
 * logs raw (t, u_cmd_L, u_cmd_R, u_act_L, u_act_R, vL, vR, iL, iR) to
 * /usd/sysid_<mode>.csv.  Then run sysid_fit.py offline to get kV, kA, kS.
 */
class SystemIdentification {
public:
  explicit SystemIdentification(Chassis &chassis);

  /**
   * @brief High-level session runner. Logs monotonic data across all phases.
   */
  void runSession(bool isAngular, const char *filename);

private:
  Chassis &chassis_;
  FILE *fp = nullptr;
  double _currentTime = 0.0;
  static constexpr double DT = 0.010;

  struct RawRow {
    double t;
    double uL, uR;
    double u_act_l, u_act_r;
    double vL, vR;
    double heading;
    const char *tag;
  };

  bool beginLog(const char *filename);
  void endLog();
  void logRow(const RawRow &r, const char *tag);
  RawRow sampleRow(double t, double uL, double uR);

  // Phase runners (Continuous time)
  void runMultiStep(double maxUL, double maxUR, double stepSize, int numSteps);
  void runPRBS(double maxUL, double maxUR, double holdTime, double duration);
  void runChirp(double maxUL, double maxUR, double minFreq, double maxFreq,
                double duration);
  void waitForButton();
};

// Helper functions for easy access from opcontrol
void run_linear_sysid(Chassis &chassis);
void run_angular_sysid(Chassis &chassis);

#endif
