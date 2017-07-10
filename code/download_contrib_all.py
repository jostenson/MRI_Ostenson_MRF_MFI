#!/usr/bin/env python
"""
Downloads all contributed code files using the URL and destination
folder specified in the comma-separated-varaiable './download_contrib_index.csv'
file, which should have the columns:

    filename,URL,folder,SHA1

Specified paths are assumed to be relative to the detected location of
this python script, i.e., the location of download_contrib_all.py

If the SHA1 hash of the downloaded file matches, the function increments the
reported SUCCESS_COUNT otherwise it increments the reported FAIL_COUNT

If a file is already downloaded, only the SHA1 checksum is performed.

Warning messages are generated when:
       * file already exists locally
       * file did not download
       * file SHA1 hash value does not match expected value

   If the file storage service changes, the URLs in the index file
   will be updated to reflect the new location.
"""
import os
import urllib2
import hashlib
import zipfile

# detect location of this script
pathstr_python_script = os.path.dirname(os.path.realpath(__file__))
print(pathstr_python_script)

# full path to download_contrib_index.csv
filename_index = os.path.join(pathstr_python_script, 'download_contrib_index.csv')

# initialize counts
success_count = 0
fail_count = 0

def sha1_for_file(filename, block_size=2**27):
    sha1 = hashlib.sha1()
    f = open(filename, 'rb')
    while True:
        data = f.read(block_size)
        if not data:
            break
        sha1.update(data)
    f.close()
    return sha1.hexdigest()

def file_is_verified(destination_filename_fullpath, sha1_expected, should_print_output):
    # Reports success if file exists on file-system and SHA1 matches.
    result = False
    if os.path.isfile(destination_filename_fullpath):
        if should_print_output:
            print("%s exists" % (destination_filename_fullpath) )
            # success only if SHA1 checksum matches
            sha1_download = sha1_for_file(destination_filename_fullpath)
        if sha1_expected==sha1_download:
            if should_print_output:
                print("Success! SHA1 checksum verified for %35s : %s (expected) == %s (downloaded)" % (filename_download, sha1_expected, sha1_download) )
                result = True
        else:
            if should_print_output:
                print("INVALID SHA1 checksum for            %35s: %s (expected) != %s (downloaded)" % (filename_download, sha1_expected, sha1_download) )
                result = False
    return result


with open(filename_index, 'r') as f:
    lines = f.read().splitlines()
    f.close()
    # loop over files. Skip header row of CSV.
    for line in lines[1:]:
        print("")
        # download details
        filename_download, URL, folder_relative_path, sha1_expected = line.strip().split(',')
        destination_filename_fullpath = os.path.join(pathstr_python_script, folder_relative_path, filename_download)

        success = 0
        if file_is_verified(destination_filename_fullpath, sha1_expected, True):
            print("Skipping download of verified file: %s" % destination_filename_fullpath)
            success = 1
            success_count += 1
        else:
            # download file
            print("Downloading file %s from %s" % (filename_download, URL) )
            filehandle = urllib2.urlopen(URL)

            # Open our local file for writing
            with open(destination_filename_fullpath, "wb") as local_file:
                local_file.write(filehandle.read())
                local_file.close()

                if file_is_verified(destination_filename_fullpath, sha1_expected, True):
                    success = 1
                    success_count += 1

                else:
                    fail_count += 1
                    print("%s did not successfully download" % (filename_download) )

	if success:
		# extract if zip file
    		if zipfile.is_zipfile(destination_filename_fullpath):
			print("Extracting zip file %s" % destination_filename_fullpath)
			z = zipfile.ZipFile(destination_filename_fullpath, "r")
			z.extractall(os.path.dirname(destination_filename_fullpath))
#          os.remove(destination_filename_fullpath)
		else:
			print("%s is not a zip file" % destination_filename_fullpath)

print("\nSuccessfully downloaded and verified %.2f%% : success_count = %d, fail_count = %d\n" % (100.0 * success_count/float(success_count+fail_count+0.000001), success_count, fail_count) )
