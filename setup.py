from setuptools import find_packages, setup
from typing import List

def get_requirements() -> List[str]:
    """
    This function will return list of reuqirements
    """
    requirements_list: List[str] = []
    try:
        with open("requirements.txt", "r") as file:
            # Read lines from the file
            lines = file.readlines()
            # process lines
            for line in lines:
                requirement = line.strip()
                if requirement and requirement != "-e .":
                    requirements_list.append(requirement)
    except FileNotFoundError:
        print("requirements.txt not found")
    return requirements_list


setup(
    name="BioMarkr",
    version = "0.0.1",
    author="Venkatesh Munaga",
    author_email="munagavenkatesh02@gmail.com",
    packages=find_packages(),
    install_required=get_requirements()
)