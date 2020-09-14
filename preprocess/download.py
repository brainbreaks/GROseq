import tarfile
import urllib
import requests
import argparse
import zipfile
import tqdm
import tempfile
import re
import shutil

def download_file(url):
	if re.match("^(http|ftp)", url):
		response = requests.get(url, stream=True)
		total_size_in_bytes= int(response.headers.get('content-length', 0))
		block_size = 1024 #1 Kibibyte

		progress_bar = tqdm.tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
		with tempfile.NamedTemporaryFile(mode="wb", delete=False) as file:
			path = file.name
			for data in response.iter_content(block_size):
				progress_bar.update(len(data))
				file.write(data)
			progress_bar.close()
		if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
			print("ERROR, something went wrong")
	else: # local file
		tmp = tempfile.NamedTemporaryFile(delete=False)
		tmp.close()
		shutil.copy2(url, tmp.name)
		path = tmp.name

	return path

def download_raw_genome(url):
	download_file(url)
	stream = urllib.request.urlopen(url)
	tar = tarfile.open(fileobj=stream, mode="r|gz")
	tar.extractall()


def download_bowtie2_index(url, dest):
	print("Downloading bowtie index '{}'\n".format(url))
	path = download_file(url)

	print("Extracting bowtie2 index into '{}'\n".format(dest))
	with zipfile.ZipFile(path, 'r') as zf:
		for member in tqdm.tqdm(zf.infolist(), desc='Extracting '):
			try:
				zf.extract(member, dest)
			except zipfile.error as e:
				pass
		# zip.extractall(dest)

def main():
	# path = download_bowtie2_index("https://file-examples-com.github.io/uploads/2017/02/zip_2MB.zip")
	# download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz")
	# download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz")
	#
	# download_bowtie2_index("https://genome-idx.s3.amazonaws.com/bt/mm9.zip", "data/mm9")
	# download_bowtie2_index("https://genome-idx.s3.amazonaws.com/bt/hg19.zip", "data/mm9")
	download_bowtie2_index("/home/s215v/Downloads/mm9.zip", "data/mm9")
	download_bowtie2_index("/home/s215v/Downloads/hg19.zip", "data/hg19")


	return

if __name__ == "__main__":
    # execute only if run as a script
    main()