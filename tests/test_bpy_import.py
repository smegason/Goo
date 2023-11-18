import pytest


@pytest.mark.parametrize(
    "imported", ("bpy"), ids=("import bpy")
)
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