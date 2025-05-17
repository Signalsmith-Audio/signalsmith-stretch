from http.server import HTTPServer, SimpleHTTPRequestHandler
import ssl
import os

file_path = os.path.dirname(__file__)
certfile = os.path.join(file_path, 'key.pem')
print(f"Using certificate file: {certfile}")

# Check if the certificate file exists and is readable
if not os.path.isfile(certfile):
    raise FileNotFoundError(f"Certificate file not found: {certfile}")

print("Starting HTTPS server...")
# Create an HTTP server with SSL
httpd = HTTPServer(('localhost', 8081), SimpleHTTPRequestHandler)
httpd.socket = ssl.wrap_socket(
  httpd.socket,
  certfile=certfile,
  server_side=True
)
httpd.serve_forever()