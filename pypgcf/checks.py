"""
Perform necessary checks during environment setup
"""
from pathlib import Path

def check_if_dir_is_empty(directory: Path):
    """
    Check if a directory is empty
    """
    items = directory.glob("*")
    if not any(items):
        return True
    return False

def check_if_dir_exists(directory: Path):
    """
    Check if a directory exists
    """
    if directory.exists():
        return True
    return False
