#!/usr/bin/env python3
"""Minimal streamable-HTTP MCP client for the SSZ Physics RAG Public server.

Used by tools/evaluate_true_closure.py to pull cited provenance for every
physics identity in the TRUE FULL CLOSURE dependency graph.  Tool names are
DISCOVERED via tools/list — never invented.
"""
from __future__ import annotations

import json
import urllib.request

ENDPOINT = "https://physics.exit-matrix.net/mcp"


def _post(payload: dict, sid: str | None = None, timeout: int = 60):
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json, text/event-stream",
    }
    if sid:
        headers["mcp-session-id"] = sid
    req = urllib.request.Request(
        ENDPOINT, data=json.dumps(payload).encode(), headers=headers,
        method="POST")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read().decode()
        new_sid = resp.headers.get("mcp-session-id", sid)
    data = None
    for line in raw.splitlines():
        if line.startswith("data:"):
            data = json.loads(line[5:])
            break
    if data is None and raw.strip():
        data = json.loads(raw)
    return data, new_sid


class PhysicsRAG:
    """Authenticated-by-session client; re-initializes when a session dies."""

    def __init__(self) -> None:
        self.sid: str | None = None
        self.server_info: dict = {}
        self._init()

    def _init(self) -> None:
        data, sid = _post({
            "jsonrpc": "2.0", "id": 1, "method": "initialize",
            "params": {
                "protocolVersion": "2024-11-05",
                "capabilities": {},
                "clientInfo": {"name": "ssz-true-closure", "version": "1.0"},
            },
        })
        self.sid = sid
        self.server_info = data["result"]["serverInfo"]
        _post({"jsonrpc": "2.0", "method": "notifications/initialized"},
              sid=self.sid)

    def call(self, name: str, arguments: dict) -> dict:
        """Call a tool; returns the parsed result content (dict)."""
        payload = {
            "jsonrpc": "2.0", "id": 2, "method": "tools/call",
            "params": {"name": name, "arguments": arguments},
        }
        try:
            data, _ = _post(payload, sid=self.sid)
        except Exception:
            self._init()
            data, _ = _post(payload, sid=self.sid)
        if "error" in data:
            raise RuntimeError(f"MCP error: {data['error']}")
        result = data["result"]
        content = result.get("content", [])
        out = []
        for item in content:
            if item.get("type") == "text":
                try:
                    out.append(json.loads(item["text"]))
                except json.JSONDecodeError:
                    out.append(item["text"])
        return {"isError": result.get("isError", False), "content": out}

    def search(self, query: str, limit: int = 5) -> dict:
        return self.call("search_physics_corpus",
                         {"query": query, "limit": limit})
