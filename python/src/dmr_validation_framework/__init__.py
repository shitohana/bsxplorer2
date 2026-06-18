"""Standalone DMR validation framework namespace.

The top-level package intentionally avoids eager imports so that importing
``dmr_validation_framework`` does not pull pandas-heavy IO, model, or workflow
modules. Import concrete helpers from their implementation modules, for
example ``dmr_validation_framework.core.io``.
"""

__all__: list[str] = []
