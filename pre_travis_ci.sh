# Prepair test settings.py file.

# Change settings.py for testing.
cp settings_test.py settings.py

# Empty output file names container.
echo -n '' > names.txt

# Run by feeding test input for further testing.
python run_sample.py
