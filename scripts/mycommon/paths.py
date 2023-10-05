from pathlib import Path

from mycommon.events import split_event_label
from mycommon.label import get_event_type_label


def get_event_label_from_path(path):
    return Path(path).parent.name


def check_same_event_type(input):
    def get_event_type_label_from_path(file):
        event_label = get_event_label_from_path(file)
        event, _ = split_event_label(event_label)
        return get_event_type_label(event)

    event_types = [get_event_type_label_from_path(file) for file in input]
    return len(set(event_types)) == 1
