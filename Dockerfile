# ####################################################################
#
# Pipeline/Dockerfile 
#
# ####################################################################
# 
# Build docker container to execute the whole IHG pipeline
# contains slave and master SGE and client and server MYSQL
# How to achieve it is still a mistery! :)
#
# Last edited		2020.02.27
# Created 		2020.02.21
#

# ####################################################################


# Inherits base container from opensuse/leap
FROM opensuse/leap


# Main working directory (when logging in into the container)
WORKDIR /sequencing


# Define frequently used variables/paths

	# Custom software directories
	ARG	PACKAGES=/usr/local/packages
	ARG	PIPELINE=/pipeline
	ARG	SEQPACKAGES=/usr/local/packages

	# System directories
	ARG	BIN=/usr/bin
	ARG	INC=/usr/include
	ARG	LIB=/usr/lib
	ARG	LOCALBIN=/usr/local/bin

# Create directory structure
	RUN mkdir -p /data
	RUN mkdir -p $PIPELINE
	RUN mkdir -p $PACKAGES
	RUN mkdir -p $SEQPACKAGES

# Install required packages
RUN zypper install -y 		\
	bzip2			\
	cmake			\
	curl			\
	curl-devel		\
	docker			\
	fontconfig-devel	\
	freetype-devel		\
	gcc			\
	gcc-c++			\
	gd-devel		\
	gdlib			\
	git			\
	gnu_parallel		\
	gzip			\
	ImageMagick		\
	ImageMagick-devel	\
	libbz2-devel		\
	libopenblas_pthreads-devel \
	libX11-devel		\
	libXext-devel		\
	libXft-devel		\
	libXpm-devel		\
	mariadb-client 		\
	ncurses-devel		\
	perl-CPAN* 		\
	perl-DateTime* 		\
	perl-DBI   		\
	perl-DBD-mysql		\
	perl-List-MoreUtils	\
	perl-Log-Log4perl 	\
	perl-Log-Dispatch*	\
	perl-MailTools		\
	perl-MIME-Lite		\
	perl-Test-Memory-Cycle	\
	perl-Tie-IxHash		\
	perl-XML-DOM-XPath	\
	perl-XML-LibXML		\
	perl-XML-Simple 	\
	python			\
	python-devel		\
	python-pip		\
	R-base 			\
	tar			\
	unzip			\
	wget			\
	xz-devel		\
	zlib-devel


	#Package amendments : Incorporate above
	#
	# RUN zypper install ..



	#################################################################################################################################################################
	# Bundled Tools
	#################################################################################################################################################################

	# Install the jdks 1.7.0 and 1.8.0
	# Since JDKs are behind a login wall and newer versions have paid license we redistribute 7 and 8 versions.
	# 1.7.0_79 and 1.8.0_65 are freely redistributable
	# 	1.7.0_79 LICENSE: https://www.oracle.com/technetwork/java/javase/jdk-7-readme-429198.html#redistribution
	# 	1.8.0_65 LICENSE: https://www.oracle.com/java/technologies/javase/jdk8-readme.html#redistribution

 	ADD build/packages /build
	RUN \
		tar -xvf /build/jdk1.7.0_79.tar.gz -C $SEQPACKAGES/;													\
		tar -xvf /build/jdk1.8.0_65.tar.gz -C $SEQPACKAGES/;	




	#################################################################################################################################################################
	# Downloadable tools
	#################################################################################################################################################################

	# Install bamUtil
	# v.1.0.14 - requires libStatGen v.1.0.14
		RUN \
			curl -L 'https://github.com/statgen/bamUtil/archive/v1.0.14.tar.gz' > $SEQPACKAGES/bamUtil-1.0.14.tar.gz;					\
			tar -xvf $SEQPACKAGES/bamUtil-1.0.14.tar.gz -C $SEQPACKAGES/;											\
			curl -L 'https://github.com/statgen/libStatGen/archive/v1.0.14.tar.gz' > $SEQPACKAGES/bamUtil-1.0.14/libStatGen.tar.gz;				\
			tar -xvf $SEQPACKAGES/bamUtil-1.0.14/libStatGen.tar.gz -C $SEQPACKAGES/bamUtil-1.0.14/;								\
			mv $SEQPACKAGES/bamUtil-1.0.14/libStatGen-1.0.14 $SEQPACKAGES/bamUtil-1.0.14/libStatGen;							\
			cd $SEQPACKAGES/bamUtil-1.0.14/;														\
			sed -e 's/-Werror//' -i Makefile */Makefile */*/Makefile */*/*/Makefile */*/*/*/Makefile;							\
			make -C $SEQPACKAGES/bamUtil-1.0.14 LIB_PATH_BAM_UTIL=$SEQPACKAGES/bamUtil-1.0.14/libStatGen;							\
			ln -s $SEQPACKAGES/bamUtil-1.0.14 $SEQPACKAGES/bamUtil;												\
			rm $SEQPACKAGES/bamUtil-1.0.14.tar.gz;														\
			rm $SEQPACKAGES/bamUtil-1.0.14/libStatGen.tar.gz;

	# Install bcftools
	#req 1.6 - version 0.1.19 embedded with samtools
		RUN \
			curl -L 'https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2' > $SEQPACKAGES/bcftools-1.6.tar.bz2;			\
			tar -xvf $SEQPACKAGES/bcftools-1.6.tar.bz2 -C $SEQPACKAGES/;											\
			cd $SEQPACKAGES/bcftools-1.6;															\
			./configure;																	\
			make -C $SEQPACKAGES/bcftools-1.6;														\
			rm  $SEQPACKAGES/bcftools-1.6.tar.bz2;

	# Install BEDTools
	# v.2.23.0
		RUN \
			curl -L 'https://github.com/arq5x/bedtools2/archive/v2.23.0.tar.gz' > $SEQPACKAGES/bedtools2-2.23.0.tar.gz;					\
			tar -xvf $SEQPACKAGES/bedtools2-2.23.0.tar.gz -C $SEQPACKAGES/;											\
			make -C $SEQPACKAGES/bedtools2-2.23.0;														\
			ln -s $SEQPACKAGES/bedtools2-2.23.0 $SEQPACKAGES/BEDTools;											\
			ln -s $SEQPACKAGES/bedtools2-2.23.0/bin/bedtools $SEQPACKAGES/bedtools2-2.23.0/bedtools;							\
			rm $SEQPACKAGES/bedtools2-2.23.0.tar.gz;

	# Install BreakDancer
	# v.1.1.2_2013_03_08 - requires samtools v0.1.16
		RUN \
			curl -L 'https://github.com/zhangjutao/Breakdancer-1.1.2/blob/master/breakdancer-1.1.2_2013_03_08.zip?raw=true' > $SEQPACKAGES/breakdancer-1.1.2.zip; \
			unzip $SEQPACKAGES/breakdancer-1.1.2.zip -d $SEQPACKAGES/;											\
			curl -L 'https://master.dl.sourceforge.net/project/samtools/samtools/0.1.6/samtools-0.1.6.tar.bz2' > $SEQPACKAGES/breakdancer-1.1.2/samtools-0.1.6.tar.bz2; \
			tar -xvf $SEQPACKAGES/breakdancer-1.1.2/samtools-0.1.6.tar.bz2 -C $SEQPACKAGES/breakdancer-1.1.2/;						\
			echo -e "all:\n\tg++ -g -Wall -O2 -I$SEQPACKAGES/breakdancer-1.1.2/samtools-0.1.6 BreakDancerMax.cpp AlnParser.cpp Poisson.cpp -o breakdancer-max -lm -lz -L$SEQPACKAGES/breakdancer-1.1.2/samtools-0.1.6 -lbam" > $SEQPACKAGES/breakdancer-1.1.2/cpp/Makefile; \
		        sed -i 's/lcurses/lncurses/g' $SEQPACKAGES/breakdancer-1.1.2/samtools-0.1.6/Makefile;								\
                        sed -i 's/\-g\ \-Wall\ \-O2/\-g\ \-Wall\ \-O2\ \-fPIC/g' $SEQPACKAGES/breakdancer-1.1.2/samtools-0.1.6/Makefile;				\
			make -C $SEQPACKAGES/breakdancer-1.1.2/samtools-0.1.6;												\
			make -C $SEQPACKAGES/breakdancer-1.1.2/cpp/;													\
			cpan install Statistics::Descriptive Math::CDF;													\
			PERL_MM_USE_DEFAULT=1 perl -MCPAN -e 'install GD::Graph::histogram';										\
			ln  -s $SEQPACKAGES/breakdancer-1.1.2 $SEQPACKAGES/breakdancer;											\
			rm  $SEQPACKAGES/breakdancer-1.1.2/samtools-0.1.6.tar.bz2;											\
			rm  $SEQPACKAGES/breakdancer-1.1.2.zip;

	# Install bwa
	# v.0.7.5a
		RUN \
			curl -L 'https://liquidtelecom.dl.sourceforge.net/project/bio-bwa/bwa-0.7.5a.tar.bz2' > $SEQPACKAGES/bwa-0.7.5a.tar.bz2;			\
			tar -xvf $SEQPACKAGES/bwa-0.7.5a.tar.bz2 -C $SEQPACKAGES/ ;											\
			make -C $SEQPACKAGES/bwa-0.7.5a all;														\
			cp $SEQPACKAGES/bwa-0.7.5a/bwa $LOCALBIN/;													\
			ln -s $SEQPACKAGES/bwa-0.7.5a $SEQPACKAGES/bwa;													\
			rm $SEQPACKAGES/bwa-0.7.5a.tar.bz2;

	# Install cutadapt
	# v.1.18
		RUN \
			pip install cutadapt==1.18;

	# Install EPACTS
	# v.3.3.2
	#	RUN \
	#		curl -L 'https://github.com/statgen/EPACTS/archive/v3.3.2.tar.gz' > $SEQPACKAGES/EPACTS-3.3.2.tar.gz;						\
	#		tar -xvf $SEQPACKAGES/EPACTS-3.3.2.tar.gz -C $SEQPACKAGES;											\
	#		cd $SEQPACKAGES/EPACTS-3.3.2;															\
	#		./configure; make; make install;														\
	#		ln -s $SEQPACKAGES/EPACTS-3.3.2 $SEQPACKAGES/EPACTS;												\
	#		rm $SEQPACKAGES/EPACTS-3.3.2.tar.gz;

	# Install FastQC
	# v.0.10.1
	# ZIP extracts to FastQC named folder:
		RUN \
			curl -L 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip' > $SEQPACKAGES/fastqc_v0.10.1.zip;			\
			unzip $SEQPACKAGES/fastqc_v0.10.1.zip -d $SEQPACKAGES;												\
			rm $SEQPACKAGES/fastqc_v0.10.1.zip; 

	# Install GATK (3 and 4)
	# v.3.8.0
		RUN \
			curl -L 'https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2' > $SEQPACKAGES/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2 ;\
			tar -xvf $SEQPACKAGES/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2  -C $SEQPACKAGES/ ;								\
			mv $SEQPACKAGES/GenomeAnalysisTK-3.8-0-ge9d806836 $SEQPACKAGES/gatk-3.8.0;									\
			ln -s $SEQPACKAGES/gatk-3.8.0 $SEQPACKAGES/gatk3;												\
			rm $SEQPACKAGES/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2;
	# v.4.1.7.0
		RUN \
			curl -L 'https://github.com/broadinstitute/gatk/releases/download/4.1.7.0/gatk-4.1.7.0.zip' > $SEQPACKAGES/gatk-4.1.7.0.zip;			\
			unzip $SEQPACKAGES/gatk-4.1.7.0.zip -d $SEQPACKAGES;												\
			ln -s $SEQPACKAGES/gatk-4.1.7.0 $SEQPACKAGES/gatk;												\
			ln -s $SEQPACKAGES/gatk-4.1.7.0 $SEQPACKAGES/gatk4;												\
			rm $SEQPACKAGES/gatk-4.1.7.0.zip;

	# Install HTSeq - Count
	# v.0.6.0
		RUN \
			pip install numpy;																\
			pip install htseq==0.6.0;
	
	# Install IGVTools
	# 2.3.68
		RUN \
			curl -L 'https://data.broadinstitute.org/igv/projects/downloads/2.3/igvtools_2.3.68.zip' > $SEQPACKAGES/igvtools_2.3.68.zip;			\
			unzip $SEQPACKAGES/igvtools_2.3.68.zip -d $SEQPACKAGES;												\
			rm $SEQPACKAGES/igvtools_2.3.68.zip;

	# Install lumpy
	# v.0.2.11
		RUN \
			curl -L 'https://github.com/arq5x/lumpy-sv/releases/download/0.2.11/lumpy-sv-0.2.11.tar.gz' > $SEQPACKAGES/lumpy-sv-0.2.11.tar.gz;		\
			tar -xvf $SEQPACKAGES/lumpy-sv-0.2.11.tar.gz -C $SEQPACKAGES/;											\
			make -C $SEQPACKAGES/lumpy-sv;															\
			rm $SEQPACKAGES/lumpy-sv-0.2.11.tar.gz
	
	# MACS2
	# v.2.1.2 (pip ava 2.1.2.1)
		RUN \
			pip install numpy;																\
			pip install macs2==2.1.2.1;
			

	# Install manta (illumina)
	# v0.29.6 (tag name)
		RUN \
			mkdir $SEQPACKAGES/manta; 															\
			cd $SEQPACKAGES/manta;																\
			git clone https://github.com/Illumina/manta.git	. ;												\
			mkdir build;																	\
			cd build;																	\
			cmake ..;																	\
			./configure;																	\
			make -C $SEQPACKAGES/manta/build; 														\
			make install;

	# Install picard
	# v 1.139
		RUN \
			curl -L 'https://github.com/broadinstitute/picard/releases/download/1.139/picard-tools-1.139.zip' > $SEQPACKAGES/picard-tools-1.139.zip;	\
			unzip $SEQPACKAGES/picard-tools-1.139.zip -d $SEQPACKAGES;											\
			ln -s  $SEQPACKAGES/picard-tools-1.139  $SEQPACKAGES/picard;											\
			rm $SEQPACKAGES/picard-tools-1.139.zip;
			
	# Install pindel
	# v.0.2.5a7
		RUN \
			curl -L 'https://github.com/genome/pindel/archive/v0.2.5a7.tar.gz' > $SEQPACKAGES/pindel-0.2.5a7.tar.gz;					\
			tar -xvf $SEQPACKAGES/pindel-0.2.5a7.tar.gz -C $SEQPACKAGES;											\
			curl -L 'https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2' > $SEQPACKAGES/pindel-0.2.5a7/samtools-1.10.tar.bz2;\
			tar -xvf $SEQPACKAGES/pindel-0.2.5a7/samtools-1.10.tar.bz2 -C $SEQPACKAGES/pindel-0.2.5a7;							\
			cd $SEQPACKAGES/pindel-0.2.5a7/samtools-1.10;													\
			./configure;																	\
			make -j 48;																	\
			cp -r $SEQPACKAGES/pindel-0.2.5a7/samtools-1.10/htslib-1.10/htslib  $SEQPACKAGES/pindel-0.2.5a7/samtools-1.10;					\
			ls; #TODO CONTINUE																\
			sed -i 's/abs(first.PosA - first.PosA1)/abs((int)first.PosA - (int)first.PosA1)/g'     */*cpp;							\
			sed -i 's/abs(first.PosB - first.PosB1)/abs((int)first.PosB - (int)first.PosB1)/g'     */*cpp;							\
			sed -i 's/abs(second.PosA - second.PosA1)/abs((int)second.PosA - (int)second.PosA1)/g' */*cpp;							\
			sed -i 's/abs(second.PosB - second.PosB1)/abs((int)second.PosB - (int)second.PosB1)/g' */*cpp;							\
			sed -i 's/abs(firstPos - secondPos)/abs((int)firstPos - (int)secondPos)/g'             */*cpp;							\
			sed -i 's/(abs(All\[index_a\].FirstPos - All\[index_b\].FirstPos)/(abs((int)All\[index_a\].FirstPos - (int)All\[index_b\].FirstPos)/g' */*cpp;	\
			sed -i 's/abs(All\[index_a\].SecondPos - All\[index_b\].SecondPos)/abs((int)All\[index_a\].SecondPos - (int)All\[index_b\].SecondPos)/g' */*cpp;\
			sed -i 's/abs(Reads_RP\[first\].PosA - Reads_RP\[first\].PosB)/abs((int)Reads_RP\[first\].PosA - (int)Reads_RP\[first\].PosB)/g' */*cpp;	\
			sed -i 's/abs(AllSV4Genotyping\[SV_index\].PosA - AllSV4Genotyping\[SV_index\].PosB)/abs((int)AllSV4Genotyping\[SV_index\].PosA - (int)AllSV4Genotyping\[SV_index\].PosB)/g' */*cpp; \
			make; make;																	\
			ln -s $SEQPACKAGES/pindel-0.2.5a7/src/pindel $SEQPACKAGES/pindel-0.2.5a7/pindel;								\
			ln -s $SEQPACKAGES/pindel-0.2.5a7/src/pindel2vcf $SEQPACKAGES/pindel-0.2.5a7/pindel2vcf;							\
			ln -s $SEQPACKAGES/pindel-0.2.5a7/src/sam2pindel $SEQPACKAGES/pindel-0.2.5a7/sam2pindel;							\
			ln -s $SEQPACKAGES/pindel-0.2.5a7 $SEQPACKAGES/pindel;												\
			rm $SEQPACKAGES/pindel-0.2.5a7/samtools-1.10.tar.bz2 $SEQPACKAGES/pindel-0.2.5a7.tar.gz;

	# Install RNA-SeQC
	# v.1.1.8
		RUN \
			mkdir $SEQPACKAGES/RNA-SeQC;															\
			curl -L 'https://data.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_v1.1.8.jar' > $SEQPACKAGES/RNA-SeQC/RNA-SeQC_v1.1.8.jar;

	# Install samtools
	# v.0.1.19
		RUN \
			curl 'https://netcologne.dl.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2' > $SEQPACKAGES/samtools-0.1.19.tar.bz2;	\
			tar -xvf $SEQPACKAGES/samtools-0.1.19.tar.bz2 -C $SEQPACKAGES/ ; 										\
			sed -i 's/lcurses/lncurses/g' $SEQPACKAGES/samtools-0.1.19/Makefile;										\
			sed -i 's/\-g\ \-Wall\ \-O2/\-g\ \-Wall\ \-O2\ \-fPIC/g' $SEQPACKAGES/samtools-0.1.19/Makefile;							\
			make -C $SEQPACKAGES/samtools-0.1.19;														\
			cp  $SEQPACKAGES/samtools-0.1.19/*h $INC/;													\
			cp  $SEQPACKAGES/samtools-0.1.19/*o $LIB/;													\
			cp  $SEQPACKAGES/samtools-0.1.19/libbam.a $LIB/; 												\
			cp  $SEQPACKAGES/samtools-0.1.19/samtools $BIN/;												\
			ln  -s $SEQPACKAGES/samtools-0.1.19 $SEQPACKAGES/samtools;											\
			rm  $SEQPACKAGES/samtools-0.1.19.tar.bz2;

	# Install snpEFF
	# v.3.4e
		RUN \
			curl -L 'https://master.dl.sourceforge.net/project/snpeff/snpEff_v3_4_core.zip' > $SEQPACKAGES/SnpEff.zip;					\
			unzip $SEQPACKAGES/SnpEff.zip -d $SEQPACKAGES;													\
			rm $SEQPACKAGES/SnpEff.zip;

	# Install STAR / STAR Fusion
	# /usr/local/packages/seq/STAR-master/bin/Linux_x86_64/STAR
	# 2.4.2a
	# /usr/local/packages/seq/STAR-master/STAR-Fusion-0.1.1/STAR-Fusion (included)
 	# 0.1.1 Fusion
		RUN \
			curl -L 'https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz' > $SEQPACKAGES/STAR_2.4.2a.tar.gz;					\
			tar -xvf $SEQPACKAGES/STAR_2.4.2a.tar.gz -C $SEQPACKAGES;											\
			cpan install Set::IntervalTree;															\
			ln -s $SEQPACKAGES/STAR-STAR_2.4.2a $SEQPACKAGES/STAR-master;											\
			rm $SEQPACKAGES/STAR_2.4.2a.tar.gz;

	# Install Tabix Perl
	# v.0.2.5
		RUN \
			curl 'https://liquidtelecom.dl.sourceforge.net/project/samtools/tabix/tabix-0.2.5.tar.bz2' > $SEQPACKAGES/tabix-0.2.5.tar.bz2;			\
			tar -xvf $SEQPACKAGES/tabix-0.2.5.tar.bz2 -C $SEQPACKAGES;											\
			make -C $SEQPACKAGES/tabix-0.2.5;														\
			bash -c 'cd $SEQPACKAGES/tabix-0.2.5/perl; perl $SEQPACKAGES/tabix-0.2.5/perl/Makefile.PL';							\
			make -C  $SEQPACKAGES/tabix-0.2.5/perl/;													\
			make -C  $SEQPACKAGES/tabix-0.2.5/perl/ install;												\
			cp  $SEQPACKAGES/tabix-0.2.5/tabix $BIN/;													\
			rm  $SEQPACKAGES/tabix-0.2.5.tar.bz2;


	# Install VCFTools
	# v.0.1.13
	# 		/usr/local/packages/seq/VCFtools
		RUN \
			curl -L 'https://github.com/vcftools/vcftools/archive/v0.1.13.tar.gz' > $SEQPACKAGES/vcftools-0.1.13.tar.gz;					\
			tar -xvf $SEQPACKAGES/vcftools-0.1.13.tar.gz -C $SEQPACKAGES/;											\
			cd $SEQPACKAGES/vcftools-0.1.13;														\
			make;																		\
			ln -s  $SEQPACKAGES/vcftools-0.1.13  $SEQPACKAGES/vcftools;											\
			ln -s  $SEQPACKAGES/vcftools-0.1.13  $SEQPACKAGES/VCFtools;											\
			rm  $SEQPACKAGES/vcftools-0.1.13.tar.gz;

	# VerifyBamId
	# v.1.1.2
		RUN \
			curl -L 'https://github.com/statgen/verifyBamID/releases/download/v1.1.2/verifyBamIDLibStatGen.1.1.2.tgz' > $SEQPACKAGES/verifyBamID-1.1.2.tar.gz;\
			tar -xvf $SEQPACKAGES/verifyBamID-1.1.2.tar.gz -C $SEQPACKAGES;											\
			mv $SEQPACKAGES/verifyBamID_1.1.2 $SEQPACKAGES/verifyBamID-1.1.2;										\
			cd $SEQPACKAGES/verifyBamID-1.1.2;														\
			sed -e 's/-Werror//' -i Makefile */Makefile */*/Makefile */*/*/Makefile */*/*/*/Makefile;							\
			make -C $SEQPACKAGES/verifyBamID-1.1.2;														\
			ln -s $SEQPACKAGES/verifyBamID-1.1.2 $SEQPACKAGES/verifyBamID;											\
			rm $SEQPACKAGES/verifyBamID-1.1.2.tar.gz;

	# Other comments:
	# DeepVariant is installed at the first run from the docker Daemon

# Non sequencing related packages:

	# ROOT (root.cern.ch)
	# v.6.20.00
		RUN \
			curl 'https://root.cern/download/root_v6.20.00.source.tar.gz' > $SEQPACKAGES/root-6.20.00.tar.gz;						\
			tar -xvf $SEQPACKAGES/root-6.20.00.tar.gz -C $SEQPACKAGES;											\
			mkdir $SEQPACKAGES/root-6.20.00/obj; 														\
			cd    $SEQPACKAGES/root-6.20.00/obj;														\
			cmake -DPYTHON_EXECUTABLE=/usr/bin/python2 .. ;													\
			make -j 48;																	\
			echo export ROOTSYS=$SEQPACKAGES/root-6.20.00/obj > /etc/profile.d/ROOT.sh;									\
			echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib >> /etc/profile.d/ROOT.sh;								\
			rm $SEQPACKAGES/root-6.20.00.tar.gz; 

		# TODO move build out and remove temp src files
		# Permanently set ENV variables (to be sure)
		ENV ROOTSYS=$SEQPACKAGES/root-6.20.00/obj
		ENV LD_LIBRARY_PATH=$ROOTSYS/lib

# Packages with special dependencies:

	# Install CNVnator
	# v.0.3.2
	# DEPENDS ON:
	# ROOT 6.0 root.cern.ch
	# samtools 0.1.13 (ok)
	
		RUN \
			curl 'https://codeload.github.com/abyzovlab/CNVnator/tar.gz/v0.3.2' > $SEQPACKAGES/CNVnator-0.3.2.tar.gz;					\
			tar -xvf $SEQPACKAGES/CNVnator-0.3.2.tar.gz -C $SEQPACKAGES;											\
			ln -s $SEQPACKAGES/samtools $SEQPACKAGES/CNVnator-0.3.2/samtools;										\
			make -C $SEQPACKAGES/CNVnator-0.3.2;														\
			ln -s $SEQPACKAGES/CNVnator-0.3.2 $SEQPACKAGES/CNVnator;											\
			rm  $SEQPACKAGES/CNVnator-0.3.2.tar.gz;
			


# Custom packages: Final config
	
	RUN  mkdir /usr/local/packages/seq
	RUN  for TOOL in /usr/local/packages/*; do ln -s $TOOL $( echo $TOOL | sed 's/packages/packages\/seq/') ; done
	



# CPAN
# Software specific dependencies have been installed already

	# Prerequisites:
	RUN zypper install -y gd gd-devel libtiff-devel libtiff5 tiff perl-DBD-mysql libXpm-devel ImageMagick

	RUN 	cp $SEQPACKAGES/samtools-0.1.19/libbam.a	/usr/lib/libbam.a; \
		cp $SEQPACKAGES/samtools-0.1.19/bam.h		/usr/include/bam.h;
	RUN cpan install 		\
		diagnostics		#
	RUN cpan install	Bio::DB::Fasta		#
	ENV LD_LIBRARY_PATH=$SEQPACKAGES/samtools-0.1.19
	ENV SAMTOOLS=$SEQPACKAGES/samtools-0.1.19 
	RUN cpan install	Bio::DB::Sam		#
	RUN cpan install	File::chmod::Recursive

	
	RUN cpan -fi \
		File::chmod File::chmod::Recursive \
		ExtUtils::MakeMaker Test::Requires Module::Implementation Params::Validate \
		GD Test::Deep Test::Most Test Bio::DB::Sam Log::Dispatch Log::Log4perl \
		File::Path File::ReadBackwards File::HomeDir Module::Build List::MoreUtils \
		DateTime::Locale DateTime::TimeZone DateTime Pod::Usage Statistics::Descriptive \
		GD::Graph::histogram Statistics::R MIME::Lite Math::Random Math::Round GD::Graph::bars GD::* \
		GD::Graph GD::Graph::bars GD::Graph::histogram Archive::Zip Set::IntervalTree \
		Bio::Root::Version Bio::DB::GFF Bio::ASN1::EntrezGene inc::latest Data::Stag \
		IO::String ExtUtils::Command::MM Test::Harness Data::Swap Inline::C Set::IntSpan::Fast::XS \
		Statistics::Basic Statistics::Distribution Statistics::KernelExtension Statistics::RankCorrelation \
		forks forks::shared Template


# R
	# Prerequisites:
	RUN zypper install -y xz-devel java-10-openjdk-devel libxml2-devel
	RUN rpm -i /build/udunits2-2.2.26-2.fc28.x86_64.rpm /build/udunits2-devel-2.2.26-2.fc28.x86_64.rpm

	ENV R_PACKAGE_LIST="XML gplots ggplot2 gplots Cairo RColorBrewer pvclust gdata matrixStats scatterplot3d snow RMySQL ggraph BiasedUrn"
	RUN \
		for R_PACKAGE in $R_PACKAGE_LIST;										\
		do														\
			echo "install.packages(\"$R_PACKAGE\", repos='https://cloud.r-project.org');" | R --no-save;		\
		done
	# ExomeDepth wants a specific version
	RUN echo "install.packages(\"https://cloud.r-project.org/src/contrib/Archive/ExomeDepth/ExomeDepth_1.1.12.tar.gz\");" | R --no-save;

	# BioConductor
	ENV R_BIOCONDUCTOR_COMMAND="source(\"https://bioconductor.org/biocLite.R\"); biocLite"
	ENV R_BIOCONDUCTOR_LIST="DESeq2 DESeq DEXSeq edgeR annotate genefilter org.Mm.eg.db org.Hs.eg.db arrayQualityMetrics gage gageData GenomicFeatures topGO goseq Rsamtools hgu95av2.db DiffBind ChIPseeker clusterProfiler ChIPpeakAnno EnsDb.Hsapiens.v75 TxDb.Mmusculus.UCSC.mm10.knownGene EnsDb.Mmusculus.v79 RUVSeq pathview mygene randomForest spp"	
	RUN \
		for BIOCONDUCTOR_PACKAGE in $R_BIOCONDUCTOR_LIST; 								\
		do														\
			R_BIOCONDUCTOR_COMMAND = echo $( $R_BIOCONDUCTOR_COMMAND"; biocLite(\""$BIOCONDUCTOR_PACKAGE"\")" );	\
		done;														\
		echo "$R_BIOCONDUCTOR_COMMAND" | R --no-save;

# Python PIP

	RUN \
		pip install cython; \
		pip install PySam pandas python-dateutil pytz scipy six scikit-learn;


# Add Pipeline and prepare configuration:

	ADD build/pipeline /pipeline
	RUN rm /pipeline/*xml
	RUN ln -s /configuration/conf.initAnalysis.xml 	/pipeline/conf.initAnalysis.xml
	RUN ln -s /configuration/current.config.xml	/pipeline/current.config.xml


# Install SSH

	RUN zypper install -y openssh
	RUN sshd-gen-keys-start

# Install Grid Engine 
# v.2011.11.p1

	RUN \
		zypper install -y tcsh pam-devel ant junit javacc motif motif-devel libXmu-devel libXp-devel;

	#Prerequisite:
	# BerkeleyDB 4.4.20
	RUN \
 		mkdir /opt/berkeley-db;												\	
		curl -L 'https://download.oracle.com/berkeley-db/db-4.4.20.tar.gz' > /opt/berkeley-db/db-4.4.20.tar.gz;		\
 		tar -xvf /opt/berkeley-db/db-4.4.20.tar.gz -C /opt/berkeley-db/;						\
		cd /opt/berkeley-db/db-4.4.20/build_unix;									\
		../dist/configure --prefix=/opt/berkeley-db/ --enable-rpc;							\
		make -j 30;													\
		make install;

	#LIBOPENSSL
	RUN \
		cd /opt;													\
		wget https://github.com/openssl/openssl/archive/OpenSSL_1_0_2u.tar.gz;						\
		tar -xvf OpenSSL_1_0_2u.tar.gz;											\
		mv openssl-OpenSSL_1_0_2u/ libssl;										\
		cd libssl;  													\
		./config; make -j 30; make install;										\
		mkdir /opt/libssl/lib; cp /opt/libssl/*a /opt/libssl/lib;

	#And finally:		
	RUN \
		export PATH=/usr/local/packages/jdk1.7.0_79/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin;  \
		curl -L 'https://downloads.sourceforge.net/project/gridscheduler/GE2011.11p1/GE2011.11p1.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fgridscheduler%2Ffiles%2Flatest%2Fdownload&ts=1589114548&use_mirror=netcologne' > /opt/GE2011.11p1.tar.gz;									\
		tar -xvf /opt/GE2011.11p1.tar.gz  -C /opt;									\
		cd /opt/GE2011.11p1/source;											\
		mv /opt/libssl/* /opt/GE2011.11p1/source/3rdparty/openssl;							\
		replace "/build/berkeleydb" "/opt/berkeley-db/" -- aimk.site;							\
		replace "0.9.8h" "1.0.2u" -- aimk.site;										\
		sed -i 's/set OPENSSL_HOME = \/build\/openssl\-install/set OPENSSL_HOME = \/opt\/GE2011.11p1\/source\/3rdparty\/openssl/g' aimk.site;			\
		./aimk -only-depend;												\
		./scripts/zerodepend;												\
		./aimk depend;													\
		./aimk -spool-classic -no-java -no-qtcsh;									\
		mkdir /opt/GE;													\
		export SGE_ROOT=/opt/GE;											\
		scripts/distinst -local -noexit -allall -y linux-x64;								
	#	mkdir /sgespool; mkdir /sgespool/qmaster; mkdir /sgespool/general/; mkdir /sgespool/local;			\
		


		
                #sed -i 's/set OPENSSL_HOME = \/build\/openssl\-install/set OPENSSL_HOME = \/opt\/GE2011.11p1\/source\/3rdparty\/openssl\nset OPENSSL_LIB_DIR = \/opt\/libssl\nset OPENSSL_EXTRA_LIBS = \"-lkrb5 -lz \"/g' aimk.site;
		#sed -i 's/izpack.home=\/off_home\/gridengine\/gui-installer\/IzPack/izpack.home=\/opt\/GE2011.11p1\/source/g' /opt/GE2011.11p1/source/build.properties;	\

# Ports to Open:

	#EXPOSE 80
	#EXPOSE ###SGE###
	#EXPOSE	MYSQL_PORT_3306_TCP


# Define machine role and start runs accordingly

	# Install mariadb server
	RUN zypper install -y *mariadb*
	RUN cpan install DBD::mysql
	
	# Install postgre if needed and DBD::Pg
	
	# Extra packages to install:
	# solr
	# mysql server

# Other cleanup
	# Remove ImageMagick default policies that impair IM functioning
	RUN rm -rf /etc/ImageMagick*/policy.xml


# Cleanup tmp files
	RUN \
		rm -rf /build;


# Unit tests

# VIM to correct and test scripts
RUN zypper install -y vim


#######################################
# To do:
######################################
#
#>  update all curls tu curl -L (to follow redirections)
#>  ln -s  /usr/local/packages/jdk1.8.0_65 /opt/

