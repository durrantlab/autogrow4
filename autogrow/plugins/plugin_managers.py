"""
Central management of AutoGrow's plugin system.

This module provides centralized access and coordination for all of AutoGrow's
plugin managers. It defines a PluginManagers class that maintains instances of
all plugin managers and provides utilities for their setup and coordination.

A single global instance of PluginManagers is created and maintained to provide
consistent access to all plugin managers throughout the application.
"""
from typing import Dict, Any
from autogrow.plugins.registry_base import plugin_managers

def setup_plugin_managers(params: Dict[str, Any]):
    """Initialize all plugin managers"""
    plugin_managers.setup_plugin_managers(params)