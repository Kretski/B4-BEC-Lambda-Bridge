import pytest
def pytest_addoption(parser): parser.addoption("--runslow", action="store_true", default=False)
def pytest_configure(config): config.addinivalue_line("markers", "slow: slow SymPy derivation")
def pytest_collection_modifyitems(config, items):
    if not config.getoption("--runslow"):
        for it in items:
            if "slow" in it.keywords: it.add_marker(pytest.mark.skip(reason="need --runslow"))
