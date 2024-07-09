import pathlib
import pytest


@pytest.fixture(scope="session")
def resource_location() -> pathlib.Path:
    return pathlib.Path("tests/resources")


@pytest.fixture(scope="session")
def chain_name() -> dict:
    return {"a": "alpha", "b": "beta"}
