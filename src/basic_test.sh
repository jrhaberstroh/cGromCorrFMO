cd ..

make src/unittest -j 8
if [ $? -eq 0 ]; then
    ./src/unittest
else
    echo "FAILED TO BUILD: Check stderr for information from gnu make"
fi
