# This file shows an example of how to run decombinator on a small test input
# Ensure you have created a virtual environment and installed decombinator via "pip install decombinator"

echo "Example pipeline aligning alpha chains..."
decombinator pipeline -in ./input/TINY_1.fq -c a -br R2 -bl 42 -ol M13
echo "Example pipeline completed aligning alpha chains."

echo "Example pipeline aligning beta chains..."
decombinator pipeline -in ./input/TINY_1.fq -c b -br R2 -bl 42 -ol M13
echo "Example pipeline completed aligning beta chains."
