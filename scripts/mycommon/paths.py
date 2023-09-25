from pathlib import Path


def get_event_label_from_path(path):
    return Path(path).parent.name
