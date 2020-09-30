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


def download_file(url, dest=None, description=None, overwrite=False):
	if description is None:
		description = 'Downloading file "{}" ==> "{}". Already exists, skipping...'.format(url, dest)

	if not overwrite and dest and os.path.isfile(dest):
		print(description)
		return dest

	if re.match("^(http|ftp)", url):
		response = requests.get(url, stream=True)
		total_size_in_bytes= int(response.headers.get('content-length', 0))
		block_size = 1024 #1 Kibibyte

		progress_bar = tqdm.tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
		with tempfile.NamedTemporaryFile(mode="wb", delete=False) as file:
			path = file.name
			progress_bar.set_description(description)
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
		print('Downloading raw genome "{}" ==> "{}". Already exists, skipping...'.format(url, dest))
		return

	path = download_file(url, description='Downloading raw genome "{}"\n'.format(url))

	dest_chromFa = "{}_{}".format(dest, "chromFa")
	with tarfile.open(path, mode="r:gz") as tf:
		progress_bar = tqdm.tqdm(tf.getmembers(), desc="Extracting raw genome archive contents into '{}'\n".format(dest_chromFa))
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


def download_bowtie2_index(url, dest, overwrite=False):
	if not overwrite and dest and len(glob.glob(os.path.join(dest, "*.bt2"))) > 0:
		print('Downloading bowtie2 index "{}" ==> "{}". Already exists, skipping...'.format(url, dest))
		return

	path = download_file(url, description="Downloading bowtie index '{}'\n".format(url))

	with zipfile.ZipFile(path, 'r') as zf:
		for member in tqdm.tqdm(zf.infolist(), desc="Extracting bowtie2 index into '{}'\n".format(dest)):
			try:
				zf.extract(member, dest)
			except zipfile.error as e:
				pass


def main():

	# This works
	download_bowtie2_index("https://genome-idx.s3.amazonaws.com/bt/mm9.zip", "data/mm9")
	download_bowtie2_index("https://genome-idx.s3.amazonaws.com/bt/mm10.zip", "data/mm10")
	download_bowtie2_index("https://genome-idx.s3.amazonaws.com/bt/hg19.zip", "data/hg19")

	download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz", "data/mm9/mm9.fa")
	download_file("http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes", "data/mm9/mm9.chrom.sizes")
	download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz", "data/mm10/mm10.fa")
	download_file("http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes", "data/mm10/mm10.chrom.sizes")
	download_raw_genome("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz", "data/hg19/hg19.fa")
	download_file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes", "data/hg19/hg19.chrom.sizes")

	# Homer packages
	download_file("http://homer.ucsd.edu/homer/configureHomer.pl", "data/dependencies/homer/configureHomer.pl")
	download_file("http://homer.ucsd.edu/homer/data/software/homer.v4.11.1.zip", "data/dependencies/homer/homer.v4.11.1.zip")
	download_file("http://homer.ucsd.edu/homer/data/genomes/hg19.v6.4.zip", "data/dependencies/homer/hg19.v6.4.zip")
	download_file("http://homer.ucsd.edu/homer/data/genomes/mm9.v6.4.zip", "data/dependencies/homer/mm9.v6.4.zip")
	download_file("http://homer.ucsd.edu/homer/data/genomes/mm10.v6.4.zip", "data/dependencies/homer/mm10.v6.4.zip")
	download_file("http://homer.ucsd.edu/homer/data/organisms/human.v6.3.zip", "data/dependencies/homer/human.v6.3.zip")
	download_file("http://homer.ucsd.edu/homer/data/organisms/mouse.v6.3.zip", "data/dependencies/homer/mouse.v6.3.zip")
	#
	# # Download libraries used in the pipeline
	download_file("https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip", "data/dependencies/bowtie2-2.2.9-linux-x86_64.zip")
	download_file("https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2", "data/dependencies/samtools-1.6.tar.bz2")
	download_file("https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz", "data/dependencies/bedtools-2.29.1.tar.gz")
	download_file("https://anaconda.org/bioconda/gfold/1.1.4/download/linux-64/gfold-1.1.4-gsl2.2_2.tar.bz2", "data/dependencies/gfold-1.1.4-gsl2.2_2.tar.bz2")
	download_file("http://www.artfiles.org/gnu.org/gsl/gsl-2.2.1.tar.gz", "data/dependencies/gsl-2.2.1.tar.gz")
	download_file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig", "data/dependencies/wigToBigWig")
	download_file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig", "data/dependencies/bigWigToWig")
	download_file("https://github.com/bxlab/bx-python/archive/master.zip", "data/dependencies/bx-python.zip")
	download_file("https://netcologne.dl.sourceforge.net/project/rseqc/RSeQC-2.6.4.tar.gz", "data/dependencies/RSeQC-2.6.4.tar.gz")
	download_file("https://github.com/numpy/numpy/releases/download/v1.16.6/numpy-1.16.6.tar.gz", "data/dependencies/numpy-1.16.6.tar.gz")
	download_file("https://github.com/cython/cython/archive/0.29.19.tar.gz", "data/dependencies/cython-0.29.19.tar.gz")
	download_file("https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/python-nose/nose-1.0.0%20(1).tar.gz", "data/dependencies/nose-1.0.0.tar.gz")
	download_file("https://github.com/pysam-developers/pysam/archive/v0.16.0.tar.gz", "data/dependencies/pysam-0.16.0.tar.gz")





	if os.path.exists("data/homer2"):
		print("HOMER is already downloaded")
	else:
		print("Download HOMER...")
		os.system("mkdir -p data/homer2; cd data/homer2; wget http://homer.ucsd.edu/homer/configureHomer.pl; perl configureHomer.pl -install homer mm9 mm10 hg19; cd cpp; make clean ")

	return

if __name__ == "__main__":
    main()