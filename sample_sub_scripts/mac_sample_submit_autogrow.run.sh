rm -r mac_sample_submit_autogrow.output_dir/*
rm mac_sample_submit_autogrow.json.out

python RunAutogrow.py --json mac_sample_submit_autogrow.json &> mac_sample_submit_autogrow.json.out
