import tarfile
import urllib
import requests
import argparse
import zipfile
import tqdm
import tempfile
import re
import shutil
import os
import glob

def download_file(url, dest=None, overwrite=False):
	if not overwrite and dest and os.path.isfile(dest):
		print('Downloading file "{}". Already exists, skipping...'.format(dest))
		return dest

	if re.match("^(http|ftp)", url):
		response = requests.get(url, stream=True)
		total_size_in_bytes= int(response.headers.get('content-length', 0))
		block_size = 1024 #1 Kibibyte

		progress_bar = tqdm.tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
		with tempfile.NamedTemporaryFile(mode="wb", delete=False) as file:
			path = file.name
			progress_bar.set_description('Downloading file "{}"...'.format(dest))
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

	if dest:
		shutil.copy2(path, dest)
		path = dest

	return path



def download_raw_genome(url, dest, overwrite=False):
	if not overwrite and dest and os.path.isfile(dest):
		print('Downloading bowtie index "{}". Already exists, skipping...'.format(dest))
		return

	print('Downloading bowtie index "{}"\n'.format(url))
	path = download_file(url)

	dest_chromFa = "{}_{}".format(dest, "chromFa")
	with tarfile.open(path, mode="r:gz") as tf:
		progress_bar = tqdm.tqdm(tf.getmembers(), desc="Extracting bowtie2 index into '{}'\n".format(dest_chromFa))
		for member in progress_bar:
			try:
				tf.extract(member, dest_chromFa)
			except tarfile.TarError as e:
				pass

	print("Joining chromosomes into single fasta file '{}'\n".format(dest))
	with open(dest, 'w') as outfile:
		# Iterate through list
		for names in glob.glob(os.path.join(dest_chromFa, "*")):
			with open(names) as infile:
				outfile.write(infile.read())
			outfile.write("\n")

	shutil.rmtree(dest_chromFa)


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


def main():
	# This works
	# download_bowtie2_index("https://genome-idx.s3.amazonaws.com/bt/mm9.zip", "data/mm9")
	# download_bowtie2_index("https://genome-idx.s3.amazonaws.com/bt/hg19.zip", "data/hg19")

	# download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz", "data/mm9/mm9.fa")
	# download_file("http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes", "data/mm9/mm9.chrom.sizes")
	# download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz", "data/hg19/hg19.fa")
	# download_file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes", "data/hg19/hg19.chrom.sizes")
	# download_file("http://homer.ucsd.edu/homer/data/genomes/hg19.v6.4.zip", "data/homer/hg19.v6.4.zip")
	# download_file("http://homer.ucsd.edu/homer/data/genomes/mm9.v6.4.zip", "data/homer/mm9.v6.4.zip")
	# download_file("http://homer.ucsd.edu/homer/data/organisms/human.v6.3.zip", "data/homer/human.v6.3.zip")
	# download_file("http://homer.ucsd.edu/homer/data/organisms/mouse.v6.3.zip", "data/homer/mouse.v6.3.zip")

	# download_file("https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip", "data/dependencies/bowtie2-2.2.9-linux-x86_64.zip")
	# download_file("https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2", "data/dependencies/samtools-1.6.tar.bz2")
	# download_file("https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz", "data/dependencies/bedtools-2.29.1.tar.gz")
	# download_file("https://anaconda.org/bioconda/gfold/1.1.4/download/linux-64/gfold-1.1.4-gsl2.2_2.tar.bz2", "data/dependencies/gfold-1.1.4-gsl2.2_2.tar.bz2")
	# download_file("http://www.artfiles.org/gnu.org/gsl/gsl-2.2.1.tar.gz", "data/dependencies/gsl-2.2.1.tar.gz")
	# download_file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig", "data/dependencies/wigToBigWig")
	# download_file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig", "data/dependencies/bigWigToWig")
	#
	#
	# Test download raw genome
	# download_raw_genome("/home/s215v/Workspace/GROseq/data/mm9/mm9_chromFa.tar.gz", "data/mm9/mm9.fa")
	#
	download_file("https://github.com/bxlab/bx-python/archive/master.zip", "data/dependencies/bx-python.zip")

#


	return

if __name__ == "__main__":
    main()