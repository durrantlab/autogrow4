#sphinx-apidoc -o source/ ../autogrow/

# Use these options
sphinx-apidoc -o source/ ../autogrow/ \
  -d 4 \
  -e \
  -f \
  -M \
  --private \
  -H "AutoGrow API Reference" \
  -A "Jacob Durrant" \
  -V "4.0.3"


make clean  # Clear old build
make html

