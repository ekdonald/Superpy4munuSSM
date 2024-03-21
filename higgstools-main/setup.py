from skbuild import setup
from setuptools import find_packages
import re
from packaging import version
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import get_cmake_version

with open("CMakeLists.txt", "r") as cmakelists:
    content = cmakelists.read()
    try:
        cmakeVersion = re.findall(
            r"cmake_minimum_required\(\W*VERSION\W*([0-9\.]+)\W*\)",
            content,
            re.MULTILINE | re.IGNORECASE,
        )[0]
        projectVersion = re.findall(
            r"project\(\W*HiggsTools\W*VERSION\W*([0-9\.]+)",
            content,
            re.MULTILINE | re.IGNORECASE,
        )[0]
    except IndexError:
        raise SKBuildError(
            "Could not read cmake and project version from CMakeLists.txt"
        )

# Add CMake as a build requirement if cmake is not installed or is too low a version
setup_requires = ["ninja"]
try:
    if version.parse(get_cmake_version()) < version.parse(cmakeVersion):
        setup_requires.append("cmake")
except SKBuildError:
    setup_requires.append("cmake")

setup(
    name="HiggsTools",
    version=projectVersion,
    description="HiggsPredictions, HiggsBounds, HiggsSignals --- comparing extended Higgs sectors to collider results.",
    author="The HiggsBounds Collaboration",
    license="GPLv3+",
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    cmake_install_dir="python/Higgs",
    cmake_args=["-DHiggsTools_BUILD_TESTING=OFF", "-DHiggsTools_BUILD_EXAMPLES=OFF"],
    python_requires=">=3.6",
    install_requires=["numpy", "scipy", "pandas", "requests", "matplotlib"],
    setup_requires=setup_requires,
)
