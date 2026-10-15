#!/bin/bash
# Double-click to serve this folder and open the viewer.
cd "$(dirname "$0")"
PORT=8731
echo "Serving $(pwd) on http://localhost:$PORT/"
echo "Keep this window open while you use the viewer; close it to stop the server."
python3 -m http.server $PORT >/dev/null 2>&1 &
SERVER=$!
sleep 1
open "http://localhost:$PORT/"
trap "kill $SERVER 2>/dev/null" EXIT
wait $SERVER
