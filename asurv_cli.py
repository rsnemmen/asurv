#!/usr/bin/env python3
"""Backward-compatibility shim. The real implementation is in src/asurv/cli.py."""
import sys
from asurv.cli import main
sys.exit(main())
