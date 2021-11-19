from pathlib import Path
import json

def pytest_configure(config):
    pwd = Path('.')
    config.cases = dict()
    for cases_file in pwd.glob('**/*_cases.json'):
        print(f'Adding test cases from {cases_file}')
        case_type = cases_file.name[:-11]
        with open(cases_file) as f:
            config.cases[case_type] = json.load(f)

def pytest_generate_tests(metafunc):
    for fixturename in metafunc.fixturenames:
        if fixturename.endswith('_case'):
            case_type = fixturename[:-5]
            cases = metafunc.config.cases[case_type]
            metafunc.parametrize(
                fixturename,
                cases.values(),
                ids = cases.keys()
            )
