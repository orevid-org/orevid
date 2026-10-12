"""Local server for the website (web/public) with no caching: the browser always gets the latest version
(with plain `http.server` it kept an old index.html that loaded an old app.js).

    py -3 web/serve.py [port]
"""
import functools
import http.server
import sys
from pathlib import Path

PUBLIC = Path(__file__).resolve().parent / "public"


class NoCache(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def log_message(self, *args):
        pass


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8796
    server = http.server.ThreadingHTTPServer(("127.0.0.1", port), functools.partial(NoCache, directory=str(PUBLIC)))
    print(f"http://localhost:{port}", flush=True)
    server.serve_forever()
