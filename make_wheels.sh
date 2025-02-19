#!/bin/bash
set -e -u -x

function cleanup_build {
    # Remove build directories and artifacts
    rm -rf build src/*.egg-info src/ampal.egg-info

    # Remove Cython-generated C files
    find -name "src/*.c" -delete

    # Remove other generated files (e.g., *.so)
    find -name "src/*.so" -delete

    echo "Cleaned build directory and Cython-generated files."
}

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "manylinux_2_34_x86_64" -w /io/wheelhouse/
    fi
}


# Install a system package required by our library
yum install -y findutils atlas-devel
rm -rf wheelhouse/
cleanup_build

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    # Skip the cp36-cp36m version
    if [[ "$PYBIN" == *"cp36-cp36m"* ]]; then
        echo "Skipping $PYBIN"
        continue
    fi
    # Skip the cp313-cp313t version
    if [[ "$PYBIN" == *"cp313-cp313t"* ]]; then
        echo "Skipping $PYBIN"
        continue
    fi

    cleanup_build
    "${PYBIN}/pip" install setuptools setuptools-scm cython networkx numpy hypothesis pytest
    "${PYBIN}/pip" wheel /io/ --no-deps -w wheelhouse/
done
cleanup_build

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    repair_wheel "$whl"
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    # Skip the cp36-cp36m version
    if [[ "$PYBIN" == *"cp36-cp36m"* ]]; then
        echo "Skipping $PYBIN"
        continue
    fi
    # Skip the cp313-cp313t version
    if [[ "$PYBIN" == *"cp313-cp313t"* ]]; then
        echo "Skipping $PYBIN"
        continue
    fi

    "${PYBIN}/pip" install ampal --no-index -f /io/wheelhouse
    "${PYBIN}/pytest" /io/tests
done

rm /io/wheelhouse/*-linux_x86_64.whl

# Create source distribution using Python 3.12
PYBIN="/opt/python/cp312-cp312/bin" 
"${PYBIN}/python" setup.py sdist -d /io/wheelhouse/
cleanup_build
