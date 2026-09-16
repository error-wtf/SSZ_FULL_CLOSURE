import json
from pathlib import Path

from ..types import ActionMember


def load_action_member(path: Path) -> ActionMember:
    obj = json.loads(path.read_text())
    if obj["release"] != "2026-09-16" or obj["action"]["epsilon_Y"] != 0.01:
        raise ValueError("Not the frozen production action")
    if obj["background"]["vector_branch"] != "A0prime = 0":
        raise ValueError("Not the zero-vector production branch")
    return ActionMember(obj["release"], obj["action"]["epsilon_Y"], "A0prime=0", path)
