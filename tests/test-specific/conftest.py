import json

def pytest_configure(config):
    with open('test-specific/alignment_cases.json') as f:
        config.alignment_cases = json.load(f)

def pytest_generate_tests(metafunc):
    if 'alignment_case' in metafunc.fixturenames:
        metafunc.parametrize(
            'alignment_case',
            metafunc.config.alignment_cases,
            ids = lambda case: f'{case["anchor"]}|{case["other"]}')
