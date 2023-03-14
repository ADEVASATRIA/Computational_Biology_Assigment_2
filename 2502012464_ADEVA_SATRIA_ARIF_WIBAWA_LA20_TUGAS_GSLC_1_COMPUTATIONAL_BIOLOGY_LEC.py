{
  "nbformat": 4,
  "nbformat_minor": 0,
  "metadata": {
    "colab": {
      "provenance": [],
      "authorship_tag": "ABX9TyM4OX5fKlb3qFx3Mk/EtPuC",
      "include_colab_link": true
    },
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3"
    },
    "language_info": {
      "name": "python"
    }
  },
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "view-in-github",
        "colab_type": "text"
      },
      "source": [
        "<a href=\"https://colab.research.google.com/github/ADEVASATRIA/Computational_Biology_Assigment_2/blob/master/2502012464_ADEVA_SATRIA_ARIF_WIBAWA_LA20_TUGAS_GSLC_1_COMPUTATIONAL_BIOLOGY_LEC.ipynb\" target=\"_parent\"><img src=\"https://colab.research.google.com/assets/colab-badge.svg\" alt=\"Open In Colab\"/></a>"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": 7,
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "skFZROaYqb_r",
        "outputId": "e4d2cfc1-7c7f-49d0-c812-3ea88f0184ba"
      },
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "Looking in indexes: https://pypi.org/simple, https://us-python.pkg.dev/colab-wheels/public/simple/\n",
            "Requirement already satisfied: biopython in /usr/local/lib/python3.9/dist-packages (1.81)\n",
            "Requirement already satisfied: numpy in /usr/local/lib/python3.9/dist-packages (from biopython) (1.22.4)\n"
          ]
        }
      ],
      "source": [
        "try:\n",
        "    import google.colab\n",
        "    # Running on Google Colab, so install Biopython first\n",
        "    !pip install biopython\n",
        "except ImportError:\n",
        "    pass"
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "import os\n",
        "import sys\n",
        "\n",
        "from urllib.request import urlretrieve\n",
        "\n",
        "import Bio\n",
        "from Bio import SeqIO, SearchIO, Entrez\n",
        "from Bio.Seq import Seq\n",
        "from Bio.SeqUtils import GC\n",
        "from Bio.Blast import NCBIWWW\n",
        "from Bio.Data import CodonTable\n",
        "\n",
        "print(\"Python version:\", sys.version_info)\n",
        "print(\"Biopython version:\", Bio.__version__)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "D82rFm2QraNe",
        "outputId": "3c61a2b2-67da-4a39-ab18-d815afa379ba"
      },
      "execution_count": 8,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "Python version: sys.version_info(major=3, minor=9, micro=16, releaselevel='final', serial=0)\n",
            "Biopython version: 1.81\n"
          ]
        }
      ]
    },
    {
      "cell_type": "markdown",
      "source": [
        "PEMBUATAN 2 **SEQUENCE**"
      ],
      "metadata": {
        "id": "54Wln32zrf0Y"
      }
    },
    {
      "cell_type": "code",
      "source": [
        "generic_dna_1 = 'ATGATCTCGTAA'\n",
        "generic_dna_2 = 'ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA'\n",
        "\n",
        "print(\"Sequence 1 = \", generic_dna_1)\n",
        "print(\"Sequence 2 = \", generic_dna_2)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "J3-xvlWHrjVx",
        "outputId": "1bae7400-4644-4f30-97fd-3c1df125c672"
      },
      "execution_count": 9,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "Sequence 1 =  ATGATCTCGTAA\n",
            "Sequence 2 =  ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA\n"
          ]
        }
      ]
    },
    {
      "cell_type": "markdown",
      "source": [
        "Menghitung Frekuensi Setiap Basa/Nukotida pada kedua Sequence"
      ],
      "metadata": {
        "id": "wlHFwqiCsAdd"
      }
    },
    {
      "cell_type": "code",
      "source": [
        "generic_dna_1 = 'ATGATCTCGTAA'\n",
        "generic_dna_2 = 'ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA'\n",
        "\n",
        "#Hiitung frekuensi setiap basa pada seq 1\n",
        "a_count_1 = generic_dna_1.count('A')\n",
        "t_count_1 = generic_dna_1.count('T')\n",
        "c_count_1 = generic_dna_1.count('C')\n",
        "g_count_1 = generic_dna_1.count('G')\n",
        "\n",
        "#Hitung frekuensi setiap basa pada seq 2\n",
        "a_count_2 = generic_dna_2.count('A')\n",
        "t_count_2 = generic_dna_2.count('T')\n",
        "c_count_2 = generic_dna_2.count('C')\n",
        "g_count_2 = generic_dna_2.count('G')\n",
        "\n",
        "#Cetak Hasil \n",
        "print(\"Berikut Hasil Frequency dari sequence 1\")\n",
        "print(\"Frequency of A in seq 1 = \", a_count_1)\n",
        "print(\"Frequency of T in seq 1 = \", t_count_1)\n",
        "print(\"Frequency of C in seq 1 = \", c_count_1)\n",
        "print(\"Frequency of G in seq 1 = \", g_count_1)\n",
        "\n",
        "print(\"\\n\")\n",
        "print(\"Berikut Hasil Frequency dari sequence 2\")\n",
        "print(\"Frequency of A in seq 1 = \", a_count_2)\n",
        "print(\"Frequency of T in seq 1 = \", t_count_2)\n",
        "print(\"Frequency of C in seq 1 = \", c_count_2)\n",
        "print(\"Frequency of G in seq 1 = \", g_count_2)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "WkGAFZOzsM2a",
        "outputId": "49f27d07-584c-4da1-d584-493877d94706"
      },
      "execution_count": 10,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "Berikut Hasil Frequency dari sequence 1\n",
            "Frequency of A in seq 1 =  4\n",
            "Frequency of T in seq 1 =  4\n",
            "Frequency of C in seq 1 =  2\n",
            "Frequency of G in seq 1 =  2\n",
            "\n",
            "\n",
            "Berikut Hasil Frequency dari sequence 2\n",
            "Frequency of A in seq 1 =  16\n",
            "Frequency of T in seq 1 =  9\n",
            "Frequency of C in seq 1 =  10\n",
            "Frequency of G in seq 1 =  4\n"
          ]
        }
      ]
    },
    {
      "cell_type": "markdown",
      "source": [
        "MENGHITUNG KANDUNGAN  GC , KANDUNGAN AT , DAN TITIK LELEH KEDUA SEQUENCE DNA"
      ],
      "metadata": {
        "id": "i8uah6oMxZoe"
      }
    },
    {
      "cell_type": "code",
      "source": [
        "generic_dna_1 = 'ATGATCTCGTAA'\n",
        "generic_dna_2 = 'ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA'\n",
        "\n",
        "# Mengitung kandungan GC dan kandungan AT pada seq 2\n",
        "gc_count_1 = generic_dna_1.count('G') + generic_dna_1.count('C')\n",
        "at_count_1 = generic_dna_1.count('A') + generic_dna_1.count('T')\n",
        "\n",
        "# Mengitung kandugan GC dan Kandungan AT pada seq 2\n",
        "gc_count_2 = generic_dna_2.count('G') + generic_dna_2.count('C')\n",
        "at_count_2 = generic_dna_2.count('A') + generic_dna_2.count('A')\n",
        "\n",
        "# Menghitung titik leleh pada seq 1 dan seq 2 \n",
        "# dengan Tm = (wA + xT) * 2 + (yG + zC) * 4\n",
        "tm_generic_dna_1 = (at_count_1 * 2) + (gc_count_1 * 4)\n",
        "tm_generic_dna_2 = (at_count_2 * 2) + (gc_count_2 * 4)\n",
        "\n",
        "# Cetak Hasil\n",
        "print(\"GC content of seq 1 = \", gc_count_1 / len(generic_dna_1) * 100, \"%\")\n",
        "print(\"AT content of seq 1 = \", at_count_1 / len(generic_dna_1) * 100, \"%\")\n",
        "print(\"Tm of seq 1 = \", tm_generic_dna_1, \"C\")\n",
        "\n",
        "print(\"\\n\")\n",
        "\n",
        "print(\"GC content of seq 2 = \", gc_count_2 / len(generic_dna_2) * 100, \"%\")\n",
        "print(\"AT content of seq 2 = \", at_count_2 / len(generic_dna_2) * 100, \"%\")\n",
        "print(\"Tm of seq 2 = \", tm_generic_dna_2, \"C\")"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "qmt7uXg6xkzu",
        "outputId": "4f065bd5-50f6-442e-b368-8abb5185371a"
      },
      "execution_count": 11,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "GC content of seq 1 =  33.33333333333333 %\n",
            "AT content of seq 1 =  66.66666666666666 %\n",
            "Tm of seq 1 =  32 C\n",
            "\n",
            "\n",
            "GC content of seq 2 =  35.8974358974359 %\n",
            "AT content of seq 2 =  82.05128205128204 %\n",
            "Tm of seq 2 =  120 C\n"
          ]
        }
      ]
    },
    {
      "cell_type": "markdown",
      "source": [
        "Menghitung Berat molekul dan titik didihnya "
      ],
      "metadata": {
        "id": "k79_RnDr8H_N"
      }
    },
    {
      "cell_type": "code",
      "source": [
        "# Diketahui Berat Molekul yang Digunakan \n",
        "# Adenin (A) = 33.3 \n",
        "# Timin (T) = 66.7 \n",
        "# Guanin (G) = 35.9\n",
        "# Sitotsin (C) = 82.1\n",
        "generic_dna_1 = 'ATGATCTCGTAA'\n",
        "generic_dna_2 = 'ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA'\n",
        "seq_length = 51\n",
        "\n",
        "# MENGHITUNG JUMLAH A, T , G, DAN C, PADA SEQUENCE 1\n",
        "a_count_1 = generic_dna_1.count('A')\n",
        "t_count_1 = generic_dna_1.count('T')\n",
        "g_count_1 = generic_dna_1.count('G')\n",
        "c_count_1 = generic_dna_1.count('C')\n",
        "\n",
        "# MENGHITUNG JUMLAH A, T , G, DAN C, PADA SEQUENCE 2\n",
        "a_count_2 = generic_dna_2.count('A')\n",
        "t_count_2 = generic_dna_2.count('T')\n",
        "g_count_2 = generic_dna_2.count('G')\n",
        "c_count_2 = generic_dna_2.count('C')\n",
        "\n",
        "# MENGHITUNG BERAT MOLEKUL SEQUENCE 1 DAN SEQUENCE 2\n",
        "mw_generic_dna_1 = (a_count_1 * 33.3) + (t_count_1 * 66.7) + (g_count_1 * 35.9) + (c_count_1 * 82.1)\n",
        "mw_generic_dna_2 = (a_count_2* 33.3) + (t_count_2 * 66.7) + (g_count_2 * 35.9) + (c_count_2 * 82.1)\n",
        "\n",
        "# Menghitung kandungan GC dan kandungan AT pada seq1 dan seq2\n",
        "gc_count_1 = (g_count_1 + c_count_1) / seq_length * 51\n",
        "at_count_1 = (a_count_1 + t_count_1) / seq_length * 51\n",
        "gc_count_1 = (g_count_2 + c_count_1) / seq_length * 51\n",
        "at_count_1 = (a_count_2 + t_count_2) / seq_length * 51\n",
        "\n",
        "# MENGHITUNG TITIK TINDIH SEQUENCE 1 DAN SEQUENCE 2 \n",
        "td_generic_dna_1 = 81.5 + (0.41 * gc_count_1) - (675 / seq_length)\n",
        "td_generic_dna_2 = 81.5 + (0.41 * gc_count_2) - (675 / seq_length)\n",
        "\n",
        "# CETAK HASIL\n",
        "print(\"BERIKUT MERUPAKAN HASIL PENGERJAAN \") \n",
        "print(\"Berat molekul Generic DNA 1 = {:.2f} g/mol\".format(mw_generic_dna_1))\n",
        "print(\"Berat molekul Generic DNA 2 = {:.2f} g/mol\".format(mw_generic_dna_2))\n",
        "print(\"Kandungan GC Generic DNA 1 = {:.2f}%\".format(gc_count_1))\n",
        "print(\"Kandungan GC Generic DNA 2 = {:.2f}%\".format(gc_count_2))\n",
        "print(\"Kandungan AT Generic DNA 1 = {:.2f}%\".format(at_count_1))\n",
        "print(\"Kandungan AT Generic DNA 2 = {:.2f}%\".format(at_count_2))\n",
        "print(\"Titik Tindih Generic DNA 1 = {:.2f} derajat C\".format(td_generic_dna_1))\n",
        "print(\"Titik Tindih Generic DNA 2 = {:.2f} derajat C\".format(td_generic_dna_2))\n"
      ],
      "metadata": {
        "id": "NU3LcGBD8P98",
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "outputId": "19859664-a35c-404a-a4d9-dd636cf679eb"
      },
      "execution_count": 19,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "BERIKUT MERUPAKAN HASIL PENGERJAAN \n",
            "Berat molekul Generic DNA 1 = 636.00 g/mol\n",
            "Berat molekul Generic DNA 2 = 2097.70 g/mol\n",
            "Kandungan GC Generic DNA 1 = 6.00%\n",
            "Kandungan GC Generic DNA 2 = 14.00%\n",
            "Kandungan AT Generic DNA 1 = 25.00%\n",
            "Kandungan AT Generic DNA 2 = 32.00%\n",
            "Titik Tindih Generic DNA 1 = 70.72 derajat C\n",
            "Titik Tindih Generic DNA 2 = 74.00 derajat C\n"
          ]
        }
      ]
    }
  ]
}
