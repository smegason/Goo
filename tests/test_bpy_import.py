import pytest
from testing_utils import parametrize_plugin_on_off


@pytest.mark.parametrize(
    "imported", ("bpy"), ids=("import bpy")
)
@parametrize_plugin_on_off
def test_bpy_import(testing_context, imported, plugin_args, expected_exitcode):
    with testing_context(
        {
            "tests/test_blender_import.py": f"""import pytest

def test_blender_import():
    import {imported}
"""
        }
    ) as ctx:
        _, stderr, exitcode = ctx.run(plugin_args)
        assert exitcode == expected_exitcode, stderr