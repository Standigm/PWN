# Copyright 2022 Standigm Inc.. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

"""PWN: Prioritization with a Warped Network."""

import logging
from importlib import metadata

from .method import pwn

logging.basicConfig(
    format='%(asctime)s [%(levelname)s] (%(module)s:%(lineno)d) %(message)s',
    level=logging.DEBUG if __debug__ else logging.INFO,
)
__all__ = ['pwn']
__version__ = metadata.version('pwn')
