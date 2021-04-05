#!/usr/bin/env bash
set -e
unset PYTHONPATH

rnaseq="$( cd .; pwd -P )"
venv="${rnaseq}/venv"

#conda_install () {
#  echo "    Installing $3 ...";
#  conda install --channel="$1" --prefix="$2" --yes --quiet "$3" >/dev/null
#  echo "    Successfully Installed $3.";
#}
#
#echo "Checking environment ..."
#for package in conda wget git
#do
#  if command -v "${package}" >/dev/null; then
#    echo "    ${package} ... installed";
#  else
#    echo "    Installing CLIP pipeline requires ${package} but it's not installed. Aborted!"; exit 1;
#  fi
#done
#echo "Checking environment complete."
#
#echo "Set up virtual environment for clip ...";
#
#echo "    Installing Python (3.8) ...";
#conda create --prefix="${venv}" --yes --quiet python=3.8 >/dev/null;
#echo "    Successful installed Python (3.8).";
#
#for package in r-base cython
#do
#  conda_install "conda-forge" "${venv}" "${package}"
#done
#
#for package in bedtools fastqc fastq-tools samtools=1.9 star=2.4.0j perl perl-app-cpanminus
#do
#  conda_install "bioconda" "${venv}" "${package}"
#done
#
#echo "    Installing Perl packages ... "
#"${clip_venv}/bin/cpanm" Statistics::Basic --quiet >/dev/null;
#"${clip_venv}/bin/cpanm" Statistics::Distributions --quiet >/dev/null;
#"${clip_venv}/bin/cpanm" Statistics::R --quiet >/dev/null;
#echo "    Successful installed 3 Perl packages."
#conda_install "mvdbeek" "${clip_venv}" "ucsc_tools";
#
#for package in cutadapt umi_tools
#do
#  conda_install "bioconda" "${clip_venv}" "${package}"
#done

echo "Installing ..."
cp source/bam2bigwig.py "${venv}/bin/bam2bigwig.py"
cp source/count_exon_junction_reads.pl "${venv}/bin/count_exon_junction_reads.pl"
cp source/count_intron_junction_reads.pl "${venv}/bin/count_intron_junction_reads.pl"
cp source/calculate_psi_for_gencode_exons_basic_manifest.pl "${venv}/bin/calculate_psi_for_gencode_exons_basic_manifest.pl"
sed "s|VIRTUAL_ENVIRONMENT|${venv}|" source/deseq.R > "${venv}/bin/deseq.R"
sed "s|VENV_BIN|${venv}/bin|g" source/process_bam_files_list.pl > "${venv}/bin/process_bam_files_list.pl"

rnaseq_py="${rnaseq}/source/rnaseq.py"
bin_rnaseq_py="${venv}/bin/rnaseq.py"
sed "s|RNASEQ_ENVIRONMENT|${venv}/environment.sh|" "${rnaseq_py}" > "${bin_rnaseq_py}"
sed -i "s|BAM2BIGWIG_PYTHON_SCRIPT|${venv}/bin/bam2bigwig.py|" "${bin_rnaseq_py}"
sed -i "s|DESEQ_R_SCRIPT|${venv}/bin/deseq.R|" "${bin_rnaseq_py}"
sed -i "s|PROCESS_BAM_FILES_LIST_PERL_SCRIPT|${venv}/bin/process_bam_files_list.pl|" "${bin_rnaseq_py}"
sed -i "s|CALCULATE_PSI|${venv}/bin/calculate_psi_for_gencode_exons_basic_manifest.pl|" "${bin_rnaseq_py}"

read -r -d '' environment << EOF || true
#!/usr/bin/env bash

unset PYTHONPATH
export PATH="${venv}/bin:\$PATH"
EOF
echo "${environment}" > "${venv}/environment.sh"

read -r -d '' script << EOF || true
#!/usr/bin/env bash

source ${venv}/environment.sh
python ${venv}/bin/rnaseq.py \$@
EOF
echo "${script}" > "${rnaseq}/rnaseq"
chmod +x "${rnaseq}/rnaseq"
echo "Successfully finalized installation."

#echo "Successfully installed and set up environment for rnaseq pipeline."
echo "Run ${rnaseq}/rnaseq -h to see the usage."
