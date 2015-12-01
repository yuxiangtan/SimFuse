#!/bin/env bash

mkdir -p $PREFIX/bin
mkdir -p $PREFIX/lib/SimFuse

cat <<EOF >$PREFIX/bin/SimFuse
#!/bin/env bash
SCRIPT_DIR=\$( cd "\$( dirname "\${BASH_SOURCE[0]}" )" && pwd )
SF_DIR=\$( cd \$SCRIPT_DIR/../lib/SimFuse && pwd )
python \$SF_DIR/SimFuse.py \$@
EOF
chmod +x $PREFIX/bin/SimFuse
cp ./SimFuse_v1/* $PREFIX/lib/SimFuse/
