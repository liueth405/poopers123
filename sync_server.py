import http.server
import socketserver
import os
import socket

PORT = 8080
BIN_DIR = "./bin"

def get_local_ip():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        s.connect(('10.255.255.255', 1))
        IP = s.getsockname()[0]
    except Exception: 
        IP = '127.0.0.1'
    finally: 
        s.close()
    return IP

class SyncHandler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/list":
            self.send_response(200)
            self.send_header("Content-type", "text/plain")
            self.end_headers()
            files = []
            if os.path.exists(BIN_DIR):
                for f in os.listdir(BIN_DIR):
                    if os.path.isfile(os.path.join(BIN_DIR, f)):
                        files.append(f"bin/{f}")
            if os.path.exists("project.pros"):
                files.append("project.pros")
            self.wfile.write("\n".join(files).encode())
        else:
            super().do_GET()

if __name__ == "__main__":
    IP = get_local_ip()
    print(f"Sync Server Active at http://{IP}:{PORT}")
    with socketserver.TCPServer(("", PORT), SyncHandler) as httpd:
        httpd.serve_forever()
