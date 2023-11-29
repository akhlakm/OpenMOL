#!/usr/bin/env bash
setup() {
    pip install --upgrade pip
    pip install --upgrade twine build
}

publish() {
    python -m build
    echo "Please use username: __token__ and password: API TOKEN"
    python -m twine upload dist/*
}

bump() {
    ## Bump the version number
    grep __version__ openmol/__init__.py
    VERSION=$(sed -n 's/__version__ = "\(.*\)"/\1/p' openmol/__init__.py)
    VERSION=$(python -c "v='$VERSION'.split('.');print('%s.%s.%d' %(v[0], v[1], int(v[2])+1))")
    echo "   >>>"
    sed -i "s/\(__version__ = \"\)[^\"]*\"/\1$VERSION\"/" openmol/__init__.py
    grep __version__ openmol/__init__.py
    git add openmol/__init__.py
    # git commit -m "Bump to version $VERSION"
}

tag() {
    # create a new git tag using the pyproject.toml version
    # and push the tag to origin
    version=$(sed -n 's/__version__ = "\(.*\)"/\1/p' openmol/__init__.py)
    git tag v$version && git push origin v$version
}

clean() {
    /bin/rm -rf dist 
}

"$@"
