"""Project-relative configuration and paths; importing this never runs analyses."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent

def local_path(value):
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = ROOT / path
    path = path.resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError(f'Path must stay inside this project: {path}')
    return path

def load_config(path=None):
    default = ROOT / 'config/default.json'
    config = json.loads(default.read_text())
    host = ROOT / 'config/host.local.json'
    if host.exists():
        config.update(json.loads(host.read_text()))
    if path is not None and Path(path).resolve() != default:
        config.update(json.loads(Path(path).read_text()))
    return config

def input_file(folder, name):
    path = Path(folder) / name
    compressed = path.with_name(path.name + '.gz')
    return compressed if not path.exists() and compressed.exists() else path

def check_output_tree(path):
    """Reject pre-existing symlinks that would redirect generated files outside."""
    path = local_path(path)
    if path.exists():
        for child in path.rglob('*'):
            if child.is_symlink():
                local_path(child)
    return path


def authorized_input_path(value):
    """Resolve authorized inputs without copying them into the public repository."""
    path = Path(value).expanduser()
    return (ROOT / path).resolve() if not path.is_absolute() else path.resolve()
