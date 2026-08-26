"""Unit tests for eefinder.get_databases (no network / no datasets binary)."""

from __future__ import annotations

import json

import pandas as pd
import pytest
from click.testing import CliRunner

from eefinder.scripts.main import cli
from eefinder.get_databases import (
    METADATA_COLUMNS,
    TRACKING_COLUMNS,
    TaxonomyRecord,
    UNINFORMATIVE_PRODUCTS,
    _genus_family_from_lineage,
    build_download_command,
    build_metadata_frame,
    cluster_identical_proteins,
    concat_protein_faas,
    count_fasta_records,
    filter_fasta_by_ids,
    find_data_report,
    find_data_reports,
    molecule_type_for_family,
    organisms_released_after,
    parse_protein_header,
    parse_release_dates,
    parse_taxonomy_report,
    reconcile_tracking,
    standardize_protein,
    validate_date,
    write_tracking_table,
)
from eefinder.normalization import strip_bracket_tags

from conftest import binaries_available


def test_standardize_protein_maps_rdrp_synonyms_within_rna_scope():
    # RdRp synonyms and spelling/spacing/case variants collapse to "RdRp".
    for raw in [
        "RNA-dependent RNA polymerase",
        "RNA dependent RNA polymerase",
        "RNA-directed RNA polymerase",
        "RdRp",
        # Compound names that previously escaped standardisation:
        "NS5B RNA-dependent RNA polymerase",
        "RNA-dependent RNA polymerase NS5",
        "RNA-dependent RNA polymerase NS5B",
        "RNA-dependent RNA polymerase protein",
        "RNA-dependent RNA polymerase P1",
        "RNA-dependent RNA polymerase P1-P2 fusion protein",
        "RNA-dependent RNA polymerase fusion protein",
        "RNA-dependent RNA polymerase catalytic subunit",
        "RNA-directed RNA polymerase nsP4",
        "3D RdRP",
        "RdRP 2",
        "Rdrp",
        "RDRP frameshift protein",
    ]:
        assert standardize_protein(raw, "ssRNA(+)") == "RdRp", raw
    # "L" / "large protein" are RdRp only for (-)ssRNA viruses.
    assert standardize_protein("L", "ssRNA(-)") == "RdRp"
    assert standardize_protein("large protein", "ssRNA(-)") == "RdRp"
    # "L polymerase" compounds on (-)RNA viruses.
    assert standardize_protein("L polymerase", "ssRNA(-)") == "RdRp"
    assert standardize_protein("L polymerase protein", "ssRNA(-)") == "RdRp"
    assert standardize_protein("L polymerase RdRp", "ssRNA(-)") == "RdRp"
    assert standardize_protein("RNA-dependent RNA polymerase L", "ssRNA(-)") == "RdRp"


def test_standardize_protein_respects_molecule_type_scope():
    # "RNA polymerase sigma factor" (phage) must NOT become RdRp, and a DNA-virus
    # "polymerase" is out of the RNA scope, so it is only cleaned/capitalised.
    assert standardize_protein("RNA polymerase sigma factor", "dsDNA") == (
        "RNA polymerase sigma factor"
    )
    assert standardize_protein("polymerase", "dsDNA") == "Polymerase"
    # "L" outside (-)ssRNA is not RdRp (just capitalised, already upper here).
    assert standardize_protein("L", "ssRNA(+)") == "L"


def test_standardize_protein_capitalizes_leading_letter():
    # Canonical names with a lower-case first letter get capitalised.
    assert standardize_protein("phosphoprotein", "ssRNA(-)") == "Phosphoprotein"
    # An unmapped name is kept but still capitalised.
    assert standardize_protein("some novel protein", "dsDNA") == "Some novel protein"
    # Names already starting upper-case (or an abbreviation) are untouched.
    assert standardize_protein("VPg", "ssRNA(+)") == "VPg"


def test_standardize_protein_removes_special_characters_and_quotes():
    # ":,/\\?!" plus double/single quotes are stripped.
    assert standardize_protein("assembly/maturation protein", "ssRNA(+)") == (
        "Assemblymaturation protein"
    )
    assert standardize_protein("weird, name!", "dsDNA") == "Weird name"
    assert standardize_protein('"gp120"', "ssRNA-RT") == "Gp120"
    assert standardize_protein("2'-O-methyltransferase", "ssRNA(+)") == (
        "2-O-methyltransferase"
    )


def test_standardize_protein_strips_directives_and_flags_bare_ones():
    # Bare directives (any case, with/without colon) become "Unknown".
    for raw in ["CDS:", "CDS", "ORF:", "ORF", "cds:", "orf"]:
        assert standardize_protein(raw, "ssRNA(+)") == "Unknown"
    # A leading directive is stripped but the real name is kept/mapped.
    assert standardize_protein("CDS: capsid protein", "ssRNA(+)") == ("Capsid Protein")
    assert standardize_protein("ORF: some novel product", "dsDNA") == (
        "Some novel product"
    )


def test_standardize_protein_matches_rdrp_in_compound_names():
    # Every compound name that mentions RdRp collapses to "RdRp" (RNA scope).
    for raw in [
        "1B protein/RdRp",
        "CP-RdRp fusion protein",
        "P2-RdRp",
        "RdRp protein",
        "RdRp-like protein",
        "putative RdRp-complex",
        "RNA-dependent RNA polymerase (RdRp)",
        "CDS: L; RNA dependent RNA polymerase; RdRp",
        "ORF 1b; RdRp; viral polymerase",
        "CDS: RNA dependant RNA polymerase; RdRp",  # misspelling
        "RdRp P1/P2 fusion protein",
        "RdRp polyprotein replicase",
        "putative CP/RdRp fusion protein",
        "putative RNA dependent RNA polymerase; RdRp",
        "putative RNA-dependent RNA polymerase RdRp",
        "RNA dependent RNA polymerase RdRp P1-P2",
        "RNA-dependent RNA-polymerase (RdRp)",
        "CDS: RdRp; RNA-dependent RNA polymerase",
        "CDS: RdRp; RNA-dependent RNA polymerase major subunit; large protein",
        "L polymerase RdRp",
    ]:
        assert standardize_protein(raw, "ssRNA(+)") == "RdRp", raw
    # ...but only within RNA scope: a dsDNA phage name is untouched.
    assert standardize_protein("RNA polymerase sigma factor", "dsDNA") == (
        "RNA polymerase sigma factor"
    )


def test_standardize_protein_maps_nucleoprotein_compounds():
    """Compound nucleoprotein names collapse to the canonical form."""
    for raw in [
        "nucleocapsid",
        "nucleocapsid phosphoprotein",
        "nucleocapsid phosphoprotein 2",
        "nucleocapsid protein",
        "nucleocapsid protein CP17",
        "nucleocapsid protein N",
        "nucleocapsid protein P8",
        "nucleocapsid shell protein",
        "nucleoprotein",
        "nucleoprotein 1",
        "nucleoprotein-like protein",
        "nucleoprotein-like protein NP",
        "pre-histone-like nucleoprotein",
        "putative DNA-binding protein major component of a nucleoprotein com",
        "putative nucleoprotein",
        # Any "...nucleocapsid..." / "...nucleoprotein..." compound collapses.
        "50 kDa protein nucleocapsid protein",
        "N nucleocapsid protein",
        "Nucleocapsid protein p3",
        "Nucleocapsid protein; encapsidates RNA genome",
        "Nucleoprotein N",
        "Structural nucleocapsid protein",
        "Structural nucleocapsid protein N",
    ]:
        for mol in ["ssRNA(-)", "ssRNA(+)", "dsDNA", ""]:
            assert standardize_protein(raw, mol) == "Nucleocapsid Protein", raw


def test_standardize_protein_normalizes_molecular_weight():
    """Molecular weight strings like '100 kDa' or '33K-like protein' are formatted."""
    assert standardize_protein("100 kDa", "dsDNA") == "100 kDa protein"
    assert standardize_protein("100 kDa protein", "dsDNA") == "100 kDa protein"
    assert standardize_protein("100K protein", "dsDNA") == "100 kDa protein"
    assert standardize_protein("100KDa protein", "dsDNA") == "100 kDa protein"
    assert standardize_protein("33K", "dsDNA") == "33 kDa protein"
    assert standardize_protein("33K protein", "dsDNA") == "33 kDa protein"
    assert standardize_protein("33K-like protein", "dsDNA") == "33 kDa protein"
    assert (
        standardize_protein("33KD putative nonstructural protein", "dsDNA")
        == "33 kDa protein"
    )
    assert standardize_protein("33L protein", "dsDNA") == "33 kDa protein"
    assert standardize_protein("33k", "dsDNA") == "33 kDa protein"
    assert standardize_protein("33k protein", "dsDNA") == "33 kDa protein"
    assert standardize_protein("33kDa protein", "dsDNA") == "33 kDa protein"
    # Hyphenated unit ("33-kDa") collapses to the same form as "33 kDa".
    assert standardize_protein("33-kDa", "dsDNA") == "33 kDa protein"
    assert standardize_protein("33-kDa protein", "dsDNA") == "33 kDa protein"
    assert standardize_protein("33-kD protein", "dsDNA") == "33 kDa protein"
    # A number fused to a leading letter (a named protein) is left alone.
    assert standardize_protein("p33-kDa", "dsDNA") == "P33-kDa"


def test_standardize_protein_fixes_common_typos():
    """Frequent NCBI misspellings are corrected so variants converge."""
    # Prefix-style typo must not corrupt the already-correct word.
    assert standardize_protein("membrane protein", "dsDNA") == "Membrane protein"
    assert standardize_protein("membran protein", "dsDNA") == "Membrane protein"
    assert standardize_protein("transmembrane protein", "dsDNA") == (
        "Transmembrane protein"
    )
    # Truncated / mangled polymerase spellings, but plurals survive.
    assert standardize_protein("RNA polymeras", "dsDNA") == "RNA polymerase"
    assert standardize_protein("RNA polymarase", "dsDNA") == "RNA polymerase"
    assert standardize_protein("RNA polymerases", "dsDNA") == "RNA polymerases"
    # Misspelled nucleocapsid / polyprotein / capsid.
    assert standardize_protein("nucleopasid protein", "ssRNA(-)") == (
        "Nucleocapsid Protein"
    )
    assert standardize_protein("polyprotien", "ssRNA(+)") == "Polyprotein"
    # A typo of an unmapped name still converges (dsDNA is out of the capsid map
    # scope, so it is corrected but not mapped to the canonical "Capsid Protein").
    assert standardize_protein("capside protein", "dsDNA") == "Capsid protein"
    # In scope, the corrected spelling reaches the canonical map entry.
    assert standardize_protein("capside protein", "ssRNA(+)") == "Capsid Protein"


def test_standardize_protein_maps_capsid_compounds():
    """Compound capsid / coat protein names collapse to the canonical form."""
    for raw in [
        "capsid protein",
        "putative capsid protein",
        "capsid protein precursor",
        "coat protein",
        "putative coat protein",
        "major coat protein",
        "minor coat protein",
    ]:
        assert standardize_protein(raw, "ssRNA(+)") == "Capsid Protein", raw


def test_standardize_protein_maps_glycoprotein_compounds():
    """All glycoprotein and spike variants collapse to Glycoprotein (G)."""
    for raw in [
        # Basic glycoprotein names
        "glycoprotein",
        "putative glycoprotein",
        "surface glycoprotein",
        "attachment glycoprotein",
        "membrane glycoprotein",
        "secreted glycoprotein",
        # Envelope glycoprotein variants (RNA + dsDNA)
        "envelope glycoprotein",
        "Envelope glycoprotein B",
        "Envelope glycoprotein H",
        "Major outer envelope glycoprotein",
        # EEV / poxvirus glycoproteins
        "EEV glycoprotein",
        "EEV membrane glycoprotein",
        "Palmytilated EEV membrane glycoprotein",
        "IEV and EEV membrane glycoprotein",
        # Spike proteins → Glycoprotein
        "spike protein",
        "spike glycoprotein",
        "tail spike protein",
        "Spike",
        "Spike structural protein",
        "Spike surface glycoprotein",
        # Structural / non-structural / transmembrane / numbered variants
        "Non-structural glycoprotein",
        "Non-structural transmembrane glycoprotein GNS",
        "Structural glycoprotein E2",
        "Virion transmembrane glycoprotein G",
        "Envelope surface glycoprotein G1",
        "GP4 glycoprotein",
        "Glycoprotein 2b (GP2b)",
        "Glycoprotein 5a",
        "Glycoprotein G",
        "Glycoprotein N",
        "Glycoproteins",
        # Truncated spellings
        "Glycop C",
        "Glycoprot",
        # A "...; fusion; ..." annotation still resolves to Glycoprotein because
        # it is not the exact "fusion glycoprotein"/"fusion protein" F entry.
        "Surface glycoprotein; fusion; receptor binding",
    ]:
        for mol in ["ssRNA(+)", "ssRNA(-)", "dsDNA", ""]:
            result = standardize_protein(raw, mol)
            assert result == "Glycoprotein", f"{raw!r} (mol={mol!r}) -> {result!r}"
    # Compound words (fused prefix) — caught via typo correction
    assert (
        standardize_protein("EEV membrane phosphoglycoprotein", "dsDNA")
        == "Glycoprotein"
    )
    assert (
        standardize_protein("Spike proteinglycoprotein", "ssRNA(+)") == "Glycoprotein"
    )
    # Fusion (F) proteins stay SEPARATE from the generic glycoprotein.
    assert standardize_protein("fusion protein", "ssRNA(-)") == "Fusion Protein"
    assert standardize_protein("fusion glycoprotein", "ssRNA(-)") == "Fusion Protein"
    assert (
        standardize_protein("putative fusion protein", "ssRNA(-)") == "Fusion Protein"
    )


def test_standardize_protein_strips_leading_qualifiers():
    """Leading hedging qualifiers are removed from the emitted name."""
    # Unmapped names keep their meaning but drop the qualifier.
    assert standardize_protein("putative membrane protein", "dsDNA") == (
        "Membrane protein"
    )
    assert standardize_protein("Predicted protein X", "dsDNA") == "Protein X"
    assert standardize_protein("Probable helicase protein", "dsDNA") == (
        "Helicase protein"
    )
    assert standardize_protein("possible small protein", "dsDNA") == "Small protein"
    # Repeated / stacked qualifiers are all removed.
    assert standardize_protein("putative probable protein X", "dsDNA") == "Protein X"
    # "hypothetical" is NOT treated as a hedging qualifier; it is dropped instead
    # (see test_standardize_protein_drops_hypothetical).


def test_standardize_protein_drops_hypothetical():
    """Any "hypothetical ..." product (incl. misspellings) becomes Unknown."""
    for raw in [
        "hypothetical protein",
        "Hypothecial protein",
        "Hyppothetical protein",
        "Hypothetical 14 kDa protein",
        "Hypothetical 4.9 kDa protein",
        "Hypothetical P1",
        "Hypothetical peptide",
        "Hypothetical polyprotein",
        "putative hypothetical protein",
    ]:
        for mol in ["ssRNA(+)", "dsDNA", ""]:
            assert standardize_protein(raw, mol) == "Unknown", raw


def test_standardize_protein_strips_ns_designation_suffix():
    """ "NSxx protein/like protein/peptide" reduces to the bare designation."""
    assert standardize_protein("NS4B protein", "ssRNA(+)") == "NS4B"
    assert standardize_protein("NS5 protein", "ssRNA(+)") == "NS5"
    assert standardize_protein("NS5-like protein", "ssRNA(+)") == "NS5"
    assert standardize_protein("NS5A peptide", "ssRNA(+)") == "NS5A"
    assert standardize_protein("NS5B protein", "ssRNA(+)") == "NS5B"
    # A bare designation is unchanged; the explicit-RdRp name still maps to RdRp.
    assert standardize_protein("NS5", "ssRNA(+)") == "NS5"
    assert standardize_protein("RNA-dependent RNA polymerase NS5", "ssRNA(+)") == "RdRp"


def test_standardize_protein_maps_capsid_variants():
    """Any "...capsid..." name collapses to Capsid Protein (RNA scope)."""
    for raw in [
        "Capsid",
        "Capsid (core) protein",
        "Capsid protein",
        "Capsid protein C",
        "Capsid prot",
    ]:
        assert standardize_protein(raw, "ssRNA(+)") == "Capsid Protein", raw
    # "nucleocapsid" must NOT be caught by the "capsid" rule (word boundary).
    assert standardize_protein("nucleocapsid protein", "ssRNA(-)") == (
        "Nucleocapsid Protein"
    )


def test_standardize_protein_maps_polymerase_plural_and_l_polyprotein():
    """Plural "polymerases" and "L polyprotein" reach RdRp within RNA scope."""
    assert standardize_protein("RNA-dependent RNA polymerases", "ssRNA(+)") == "RdRp"
    assert standardize_protein("L polyprotein", "ssRNA(-)") == "RdRp"
    # dsDNA phage "RNA polymerases" is out of RNA scope -> plural kept, not RdRp.
    assert standardize_protein("RNA polymerases", "dsDNA") == "RNA polymerases"


def test_standardize_protein_handles_misspellings():
    """Common typos in NCBI protein names are corrected before matching."""
    # "polymrease" (transposition typo)
    assert standardize_protein("RNA-dependent RNA polymrease", "ssRNA(+)") == "RdRp"
    # "polymeras" (truncated)
    assert standardize_protein("RNA-dependent RNA polymeras", "ssRNA(+)") == "RdRp"
    # "nucleocapside" (French spelling)
    assert standardize_protein("nucleocapside", "ssRNA(-)") == "Nucleocapsid Protein"


def test_standardize_protein_dispatches_per_target():
    """The ``target`` argument selects per-target standardisation logic."""
    # Virus (default) applies the viral map; RdRp synonyms collapse.
    assert standardize_protein("RNA-dependent RNA polymerase", "ssRNA(+)") == "RdRp"
    assert (
        standardize_protein("RNA-dependent RNA polymerase", "ssRNA(+)", target="virus")
        == "RdRp"
    )
    # Bacteria/host have no viral map: the same name is only cleaned/capitalised,
    # but shared cleaning (typos, capitalisation, Unknown collapse) still applies.
    assert (
        standardize_protein("RNA-dependent RNA polymerase", "", target="bacteria")
        == "RNA-dependent RNA polymerase"
    )
    assert standardize_protein("membran protein", "", target="bacteria") == (
        "Membrane protein"
    )
    assert standardize_protein("some gene", "", target="host") == "Some gene"
    for target in ("bacteria", "host"):
        assert standardize_protein("CDS:", "", target=target) == "Unknown"


def test_standardize_protein_rejects_unknown_target():
    with pytest.raises(ValueError):
        standardize_protein("capsid protein", "ssRNA(+)", target="plant")


def test_strip_bracket_tags_removes_ncbi_metadata_tags():
    """Leaked NCBI ``[key=value]`` tags are removed from protein names."""
    assert strip_bracket_tags("nucleoprotein [organism=Rabies lyssavirus]").strip() == (
        "nucleoprotein"
    )
    assert strip_bracket_tags("[gbkey=CDS] [organism=Foo virus] capsid").split() == [
        "capsid"
    ]
    # A tag value may contain one level of nested "[...]" (a strain/isolate tag).
    assert strip_bracket_tags(
        "AC1 protein [organism=Ageratum leaf curl virus - [G52]]"
    ).strip() == ("AC1 protein")
    assert strip_bracket_tags(
        "AC1 [organism=Papaya virus [IndiaPan2008]] [isolate=Panipat [INP08]]"
    ).strip() == ("AC1")
    # The header parser recovers a clean product and the nested organism name.
    header = parse_protein_header(
        "YP_1.1 AC1 protein [organism=Maize streak virus - A[South Africa]]"
    )
    assert header.product == "AC1 protein"
    assert header.organism == "Maize streak virus - A[South Africa]"
    # standardize_protein cleans the tag out of the emitted name.
    assert (
        standardize_protein("capsid protein [organism=Some virus]", "ssRNA(+)")
        == "Capsid Protein"
    )
    assert (
        standardize_protein("nucleoprotein [gene=N]", "dsDNA", target="bacteria")
        == "Nucleoprotein"
    )
    # parse_protein_header keeps the product free of leaked tags.
    header = parse_protein_header("YP_1.1 nucleoprotein [organism=Foo] [isolate=Bar]")
    assert header.product == "nucleoprotein"
    assert header.organism == "Foo"


def test_molecule_type_for_family_from_ictv_table():
    # Values from the bundled ICTV genome-composition table.
    assert molecule_type_for_family("Rhabdoviridae") == "ssRNA(-)"
    assert molecule_type_for_family("Retroviridae") == "ssRNA-RT"
    assert molecule_type_for_family("Flaviviridae") == "ssRNA(+)"
    # Unknown / bacterial families are not covered -> empty string.
    assert molecule_type_for_family("Enterobacteriaceae") == ""
    assert molecule_type_for_family("") == ""


def test_genus_family_from_lineage_uses_ictv_suffixes():
    lineage = [
        {"name": "Viruses"},
        {"name": "Kitrinoviricota"},
        {"name": "Flasuviricetes"},
        {"name": "Amarillovirales"},
        {"name": "Flaviviridae"},  # family (-viridae)
        {"name": "Orthoflavivirus"},  # genus (single word, -virus)
        {"name": "Orthoflavivirus denguei"},  # species (two words)
    ]
    genus, family = _genus_family_from_lineage(lineage)
    assert family == "Flaviviridae"
    assert genus == "Orthoflavivirus"


def test_genus_family_from_lineage_blank_when_no_ictv_names():
    genus, family = _genus_family_from_lineage(
        [{"name": "Viruses"}, {"name": "unclassified phage"}]
    )
    assert (genus, family) == ("", "")


def test_parse_protein_header_refseq_format():
    header = parse_protein_header(">YP_009664712.1 N protein [Bas-Congo tibrovirus]")
    assert header.accession == "YP_009664712.1"
    assert header.product == "N protein"
    assert header.organism == "Bas-Congo tibrovirus"


def test_parse_protein_header_datasets_cds_format():
    # NCBI datasets CDS proteins tag the organism as [organism=...] and may carry
    # extra [key=value] groups; the product is the text before the first bracket.
    header = parse_protein_header(
        ">YP_013613119.1:1-301 P1 [polyprotein=polyprotein] "
        "[organism=Paris potyvirus 4] [isolate=YLJ]"
    )
    assert header.accession == "YP_013613119.1:1-301"  # kept whole (BLAST id)
    assert header.product == "P1"
    assert header.organism == "Paris potyvirus 4"


def test_parse_protein_header_datasets_empty_product():
    header = parse_protein_header(
        ">NC_139268.1:152-1432 CDS:  [organism=Perhabdovirus trutta] [isolate=18/203]"
    )
    assert header.accession == "NC_139268.1:152-1432"
    assert header.product == "CDS:"
    assert header.organism == "Perhabdovirus trutta"


def test_parse_protein_header_without_organism():
    header = parse_protein_header("ABC123.1 hypothetical protein")
    assert header.accession == "ABC123.1"
    assert header.product == "hypothetical protein"
    assert header.organism == ""


def test_parse_protein_header_product_with_internal_brackets():
    # Only the trailing [...] is the organism; earlier brackets stay in product.
    header = parse_protein_header(">X.1 RNA polymerase [subunit 2] [Zika virus]")
    assert header.product == "RNA polymerase [subunit 2]"
    assert header.organism == "Zika virus"


def test_parse_taxonomy_report_infers_ranks_and_reads_host(tmp_path):
    # Mirrors the real NCBI datasets schema: lineage is unranked {name, taxId},
    # host is a top-level object, and there is no molType field.
    report = tmp_path / "data_report.jsonl"
    lines = [
        {
            "accession": "NC_139268.1",
            "host": {"organismName": "Percidae"},
            "virus": {
                "organismName": "Perhabdovirus trutta",
                "lineage": [
                    {"name": "Viruses", "taxId": 10239},
                    {"name": "Rhabdoviridae", "taxId": 11270},
                    {"name": "Alpharhabdovirinae", "taxId": 2810308},
                    {"name": "Perhabdovirus", "taxId": 1962501},
                    {"name": "Perhabdovirus trutta", "taxId": 1987017},
                ],
            },
        },
        {
            "accession": "NC_2",
            "virus": {
                "organismName": "uncultured phage cr116_1",
                "lineage": [{"name": "Viruses", "taxId": 10239}],
            },
        },
    ]
    report.write_text("\n".join(json.dumps(entry) for entry in lines) + "\n")

    records = parse_taxonomy_report(str(report))

    trutta = records["Perhabdovirus trutta"]
    # Family (-viridae) and genus (single word -virus) inferred by ICTV suffix.
    assert trutta.family == "Rhabdoviridae"
    assert trutta.genus == "Perhabdovirus"
    assert trutta.host == "Percidae"
    assert trutta.mol_type == ""  # never present in the datasets report

    # A non-ICTV leaf (a phage with no ranked lineage) degrades to blanks.
    phage = records["uncultured phage cr116_1"]
    assert phage.family == ""
    assert phage.genus == ""
    assert phage.host == ""


def test_build_metadata_frame_joins_headers_and_taxonomy(tmp_path):
    fasta = tmp_path / "protein.fa"
    fasta.write_text(
        ">YP_1.1 N protein [Bas-Congo tibrovirus]\nMKKV\n"
        ">YP_2.1 polymerase [Unknown virus]\nMAAA\n"
    )
    taxonomy = {
        "Bas-Congo tibrovirus": TaxonomyRecord(
            species="Bas-Congo tibrovirus",
            genus="Tibrovirus",
            family="Rhabdoviridae",
            mol_type="ssRNA(-)",
            host="Homo sapiens",
        )
    }

    frame = build_metadata_frame(str(fasta), taxonomy)

    assert list(frame.columns) == METADATA_COLUMNS
    first = frame.iloc[0]
    assert first["Accession"] == "YP_1.1"
    assert first["Species"] == "Bas-Congo tibrovirus"
    assert first["Family"] == "Rhabdoviridae"
    assert first["Protein"] == "N protein"
    assert first["Host"] == "Homo sapiens"
    # Molecule_type comes from the ICTV table via the family, not the report.
    assert first["Molecule_type"] == "ssRNA(-)"

    # An organism absent from the taxonomy still yields a row with blanks.
    second = frame.iloc[1]
    assert second["Accession"] == "YP_2.1"
    assert second["Species"] == "Unknown virus"
    assert second["Family"] == ""
    assert second["Genus"] == ""
    assert second["Molecule_type"] == ""


def test_build_metadata_frame_standardizes_when_requested(tmp_path):
    fasta = tmp_path / "protein.fa"
    fasta.write_text(">YP_1.1 RNA-dependent RNA polymerase [Some virus]\nMKKV\n")
    taxonomy = {
        "Some virus": TaxonomyRecord(
            species="Some virus",
            genus="Somevirus",
            family="Flaviviridae",  # -> Molecule_type ssRNA(+) via ICTV table
            mol_type="",
            host="",
        )
    }

    raw = build_metadata_frame(str(fasta), taxonomy, standardize=False)
    assert raw.iloc[0]["Protein"] == "RNA-dependent RNA polymerase"

    std = build_metadata_frame(str(fasta), taxonomy, standardize=True)
    # ssRNA(+) is in RdRp scope, so the synonym collapses to the canonical name.
    assert std.iloc[0]["Protein"] == "RdRp"


def test_build_metadata_frame_standardize_drops_bare_directives(tmp_path):
    fasta = tmp_path / "protein.fa"
    fasta.write_text(
        ">A.1 capsid protein [Some virus]\nMK\n"
        ">B.1 CDS: [Some virus]\nMA\n"  # bare directive -> dropped
        ">C.1 RdRp [Some virus]\nMC\n"
    )
    taxonomy = {
        "Some virus": TaxonomyRecord(
            species="Some virus",
            genus="Somevirus",
            family="Flaviviridae",  # -> ssRNA(+)
            mol_type="",
            host="",
        )
    }

    std = build_metadata_frame(str(fasta), taxonomy, standardize=True)
    # B.1 (bare "CDS:") is dropped; A.1/C.1 kept and standardised.
    assert list(std["Accession"]) == ["A.1", "C.1"]
    assert list(std["Protein"]) == ["Capsid Protein", "RdRp"]

    # Without standardisation nothing is dropped.
    raw = build_metadata_frame(str(fasta), taxonomy, standardize=False)
    assert list(raw["Accession"]) == ["A.1", "B.1", "C.1"]


def test_filter_fasta_by_ids_drops_records_not_kept(tmp_path):
    fasta = tmp_path / "db.fa"
    fasta.write_text(">A.1 x\nMK\nMK\n>B.1 y\nMA\n>C.1 z\nMC\n")

    dropped = filter_fasta_by_ids(str(fasta), {"A.1", "C.1"})

    assert dropped == 1
    text = fasta.read_text()
    assert ">A.1" in text and ">C.1" in text
    assert ">B.1" not in text and "MA" not in text
    # multi-line sequence of a kept record is preserved
    assert text.count("MK") == 2


def test_build_download_command_virus_refseq():
    command = build_download_command(
        "virus", "Flaviviridae", "/out/virus.zip", refseq=True
    )
    assert command == (
        "datasets download virus genome taxon Flaviviridae "
        "--include protein --filename /out/virus.zip --refseq"
    )


def test_build_download_command_genome_all_sequences():
    command = build_download_command(
        "host", "Aedes aegypti", "/out/host.zip", refseq=False
    )
    # A multi-word taxon is quoted; without --refseq no source filter is added.
    assert "datasets download genome taxon 'Aedes aegypti'" in command
    assert "--assembly-source" not in command
    assert command.endswith("--filename /out/host.zip")


def test_build_download_command_rejects_unknown_dataset():
    try:
        build_download_command("plant", "Zea mays", "/out/x.zip")
    except ValueError as err:
        assert "plant" in str(err)
    else:  # pragma: no cover
        raise AssertionError("expected ValueError for an unknown dataset type")


def test_concat_protein_faas_merges_all_and_counts(tmp_path):
    data = tmp_path / "ncbi_dataset" / "data"
    (data / "ASM1").mkdir(parents=True)
    (data / "ASM2").mkdir(parents=True)
    (data / "ASM1" / "protein.faa").write_text(">a\nMK\n")
    (data / "ASM2" / "protein.faa").write_text(">b\nMA\n")
    out = tmp_path / "merged.fa"

    counts = concat_protein_faas(str(tmp_path), str(out))

    assert counts.files == 2
    assert counts.total == 2
    assert counts.written == 2
    merged = out.read_text()
    assert ">a" in merged and ">b" in merged


def test_concat_protein_faas_excludes_uninformative_products(tmp_path):
    data = tmp_path / "ncbi_dataset" / "data"
    data.mkdir(parents=True)
    (data / "protein.faa").write_text(
        ">A.1 capsid protein [organism=Some virus]\nMK\n"
        ">B.1 hypothetical protein [organism=Some virus]\nMA\n"
        ">C.1 UNCHARACTERIZED protein [organism=Some virus]\nMC\n"
        ">D.1 replicase [organism=Some virus]\nMD\n"
    )
    out = tmp_path / "merged.fa"

    counts = concat_protein_faas(
        str(tmp_path), str(out), exclude_products=UNINFORMATIVE_PRODUCTS
    )

    merged = out.read_text()
    assert counts.files == 1
    assert counts.total == 4  # four records seen
    assert counts.written == 2  # two informative kept
    assert ">A.1" in merged and ">D.1" in merged  # informative kept
    assert ">B.1" not in merged  # hypothetical dropped
    assert ">C.1" not in merged  # uncharacterized dropped (case-insensitive)
    assert "MA" not in merged and "MC" not in merged  # their sequences too


def test_concat_protein_faas_raises_when_missing(tmp_path):
    try:
        concat_protein_faas(str(tmp_path), str(tmp_path / "out.fa"))
    except FileNotFoundError:
        pass
    else:  # pragma: no cover
        raise AssertionError("expected FileNotFoundError when no protein.faa exists")


def test_find_data_report(tmp_path):
    assert find_data_report(str(tmp_path)) is None
    nested = tmp_path / "ncbi_dataset" / "data"
    nested.mkdir(parents=True)
    report = nested / "data_report.jsonl"
    report.write_text("{}\n")
    assert find_data_report(str(tmp_path)) == str(report)


def test_cli_get_databases_errors_without_datasets_binary(monkeypatch, tmp_path):
    # When the datasets CLI is absent, the command fails fast with guidance.
    monkeypatch.setattr("eefinder.scripts.main.shutil.which", lambda name: None)
    result = CliRunner().invoke(
        cli,
        [
            "get-databases",
            "virus",
            "-tx",
            "Flaviviridae",
            "-od",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 1
    assert "datasets" in result.output


def test_cli_virus_defaults_taxon_to_viruses_root(monkeypatch, tmp_path):
    # The virus subcommand without -tx falls back to 10239 (Viruses).
    captured = {}

    def fake_get_databases(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(
        "eefinder.scripts.main.shutil.which", lambda name: "/usr/bin/datasets"
    )
    monkeypatch.setattr("eefinder.scripts.main.GetDatabases", fake_get_databases)

    result = CliRunner().invoke(
        cli,
        ["get-databases", "virus", "-od", str(tmp_path / "out")],
    )
    assert result.exit_code == 0, result.output
    assert captured["dataset"] == "virus"
    assert captured["taxon"] == "10239"
    assert captured["standardize_proteins"] is True


def test_cli_bacteria_defaults_taxon_to_bacteria_root(monkeypatch, tmp_path):
    # The bacteria subcommand without -tx falls back to 2 (Bacteria) and
    # standardizes (generic cleaning; there is no bacterial name map yet).
    captured = {}

    def fake_get_databases(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(
        "eefinder.scripts.main.shutil.which", lambda name: "/usr/bin/datasets"
    )
    monkeypatch.setattr("eefinder.scripts.main.GetDatabases", fake_get_databases)

    result = CliRunner().invoke(
        cli,
        ["get-databases", "bacteria", "-od", str(tmp_path / "out")],
    )
    assert result.exit_code == 0, result.output
    assert captured["dataset"] == "bacteria"
    assert captured["taxon"] == "2"
    assert captured["standardize_proteins"] is True

    # --raw-proteins opts out of the cleaning.
    captured.clear()
    result = CliRunner().invoke(
        cli,
        ["get-databases", "bacteria", "-od", str(tmp_path / "out"), "--raw-proteins"],
    )
    assert result.exit_code == 0, result.output
    assert captured["standardize_proteins"] is False


def test_cli_requires_taxon_for_host(monkeypatch, tmp_path):
    # host has no default taxon, so omitting -tx is a usage error.
    monkeypatch.setattr(
        "eefinder.scripts.main.shutil.which", lambda name: "/usr/bin/datasets"
    )
    result = CliRunner().invoke(
        cli,
        ["get-databases", "host", "-od", str(tmp_path / "out")],
    )
    assert result.exit_code != 0
    assert "taxon" in result.output


def test_cli_get_databases_rejects_unknown_subcommand(tmp_path):
    result = CliRunner().invoke(
        cli,
        ["get-databases", "plant", "-od", str(tmp_path)],
    )
    # click rejects an unknown subcommand before any command body runs.
    assert result.exit_code != 0
    assert "plant" in result.output


def test_metadata_frame_roundtrips_to_expected_csv_schema(tmp_path):
    # The written CSV must match the header the screening command expects.
    fasta = tmp_path / "protein.fa"
    fasta.write_text(">A.1 cap [Some virus]\nMK\n")
    frame = build_metadata_frame(str(fasta), {})
    csv = tmp_path / "meta.csv"
    frame.to_csv(csv, index=False)
    header = pd.read_csv(csv).columns.tolist()
    assert header == METADATA_COLUMNS


def test_count_fasta_records_counts_headers(tmp_path):
    fasta = tmp_path / "p.fa"
    fasta.write_text(">a\nMKAA\n>b\nMK\nAA\n>c\nMM\n")
    assert count_fasta_records(str(fasta)) == 3
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    assert count_fasta_records(str(empty)) == 0


@pytest.mark.skipif(not binaries_available("cd-hit"), reason="requires cd-hit on PATH")
def test_cluster_identical_proteins_removes_exact_duplicates(tmp_path):
    """100%/100% clustering drops an exact duplicate and keeps distinct ones."""
    fasta = tmp_path / "db.fa"
    # PROT_B is byte-identical to PROT_A (a duplicate); PROT_C is distinct. A
    # diverse peptide avoids cd-hit's low-complexity edge cases.
    peptide = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
    fasta.write_text(
        f">PROT_A rep [Some virus]\n{peptide}\n"
        f">PROT_B dup [Some virus]\n{peptide}\n"
        f">PROT_C other [Other virus]\n{peptide[::-1]}\n"
    )
    removed = cluster_identical_proteins(str(fasta), threads=1)

    assert removed == 1
    assert count_fasta_records(str(fasta)) == 2
    # The cd-hit .clstr sidecar is cleaned up, not left next to the database.
    assert not (tmp_path / "db.fa.nr.clstr").exists()
    assert not (tmp_path / "db.fa.nr").exists()
    ids = {
        line[1:].split()[0]
        for line in fasta.read_text().splitlines()
        if line.startswith(">")
    }
    # The representative of the identical pair survives, plus the distinct one.
    assert "PROT_C" in ids
    assert len(ids & {"PROT_A", "PROT_B"}) == 1


# --------------------------------------------------------------------------- #
# tracking.tsv: the fate of every downloaded accession
# --------------------------------------------------------------------------- #
class TestTracking:
    """Every sequence that leaves the database must say why it left."""

    def _faa(self, tmp_path, records):
        data = tmp_path / "ncbi_dataset" / "data"
        data.mkdir(parents=True)
        (data / "protein.faa").write_text(
            "".join(f">{header}\n{seq}\n" for header, seq in records)
        )
        return str(tmp_path)

    def test_uninformative_records_are_recorded_with_their_reason(self, tmp_path):
        extract = self._faa(
            tmp_path,
            [
                ("PROT_A glycoprotein [organism=Mivirus chuvi]", "MKALL"),
                ("PROT_B hypothetical protein [organism=Mivirus chuvi]", "KKTTT"),
            ],
        )
        tracking = {}
        counts = concat_protein_faas(
            extract,
            str(tmp_path / "out.faa"),
            exclude_products=UNINFORMATIVE_PRODUCTS,
            tracking=tracking,
        )

        assert counts.total == 2 and counts.written == 1
        assert tracking["PROT_A"]["Status"] == "kept"
        assert tracking["PROT_B"]["Status"] == "removed"
        assert tracking["PROT_B"]["Reason"] == "uninformative_product"
        # The dropped one is only visible here: it never reaches the FASTA.
        assert "PROT_B" not in (tmp_path / "out.faa").read_text()

    def test_every_downloaded_accession_is_present(self, tmp_path):
        extract = self._faa(
            tmp_path,
            [(f"PROT_{i} some protein [organism=X]", "MKALL") for i in range(5)],
        )
        tracking = {}
        concat_protein_faas(extract, str(tmp_path / "out.faa"), tracking=tracking)
        assert len(tracking) == 5

    def test_a_rename_is_recorded(self, tmp_path):
        fasta = tmp_path / "proteins.faa"
        fasta.write_text(
            ">PROT_A RNA-dependent RNA polymerase [organism=Mivirus chuvi]\nMKALL\n"
        )
        tracking = {
            "PROT_A": {
                "Accession": "PROT_A",
                "Species": "Mivirus chuvi",
                "Protein_downloaded": "RNA-dependent RNA polymerase",
                "Protein_final": "",
                "Name_changed": "",
                "Status": "kept",
                "Reason": "",
                "Cluster": "",
                "Cluster_representative": "",
            }
        }
        taxonomy = {
            "Mivirus chuvi": TaxonomyRecord(
                "Mivirus chuvi", "Mivirus", "Chuviridae", "", ""
            )
        }
        build_metadata_frame(str(fasta), taxonomy, standardize=True, tracking=tracking)

        row = tracking["PROT_A"]
        assert row["Protein_final"] == "RdRp"
        assert row["Name_changed"] == "yes"
        assert row["Status"] == "kept"

    def test_an_unchanged_name_says_so(self, tmp_path):
        fasta = tmp_path / "proteins.faa"
        fasta.write_text(">PROT_A RdRp [organism=Mivirus chuvi]\nMKALL\n")
        tracking = {
            "PROT_A": {
                "Accession": "PROT_A",
                "Protein_final": "",
                "Name_changed": "",
                "Status": "kept",
                "Reason": "",
            }
        }
        build_metadata_frame(str(fasta), {}, standardize=True, tracking=tracking)
        assert tracking["PROT_A"]["Name_changed"] == "no"

    def test_records_dropped_by_standardization_are_recorded(self, tmp_path):
        fasta = tmp_path / "proteins.faa"
        fasta.write_text(">PROT_A CDS: [organism=Mivirus chuvi]\nMKALL\n")
        tracking = {
            "PROT_A": {
                "Accession": "PROT_A",
                "Protein_final": "",
                "Name_changed": "",
                "Status": "kept",
                "Reason": "",
            }
        }
        frame = build_metadata_frame(
            str(fasta), {}, standardize=True, tracking=tracking
        )

        assert frame.empty
        assert tracking["PROT_A"]["Status"] == "removed"
        assert tracking["PROT_A"]["Reason"] == "product_standardized_to_unknown"

    def test_a_record_that_vanished_without_explanation_is_caught(self, tmp_path):
        """cd-hit silently drops very short sequences; the table must not lie."""
        fasta = tmp_path / "final.faa"
        fasta.write_text(">PROT_A\nMKALLVGTSGAGKSTLLQ\n")
        tracking = {
            "PROT_A": {"Accession": "PROT_A", "Status": "kept", "Reason": ""},
            "PROT_SHORT": {"Accession": "PROT_SHORT", "Status": "kept", "Reason": ""},
            "PROT_DROPPED": {
                "Accession": "PROT_DROPPED",
                "Status": "removed",
                "Reason": "uninformative_product",
            },
        }
        assert reconcile_tracking(tracking, str(fasta)) == 1

        assert tracking["PROT_A"]["Status"] == "kept"
        assert tracking["PROT_SHORT"]["Status"] == "removed"
        assert tracking["PROT_SHORT"]["Reason"] == "absent_from_final_database"
        # An already-explained removal keeps its own reason.
        assert tracking["PROT_DROPPED"]["Reason"] == "uninformative_product"

    def test_the_table_has_the_documented_columns(self, tmp_path):
        tracking = {
            "PROT_A": {column: "" for column in TRACKING_COLUMNS},
        }
        tracking["PROT_A"]["Accession"] = "PROT_A"
        tracking["PROT_A"]["Status"] = "kept"
        out = tmp_path / "tracking.tsv"
        assert write_tracking_table(tracking, str(out)) == 1

        table = pd.read_csv(out, sep="\t")
        assert list(table.columns) == TRACKING_COLUMNS


@pytest.mark.skipif(not binaries_available("cd-hit"), reason="requires cd-hit on PATH")
class TestTrackingOfDuplicates:
    """A dropped duplicate must point at the sequence that replaced it."""

    def test_duplicates_name_their_cluster_and_representative(self, tmp_path):
        fasta = tmp_path / "proteins.faa"
        # PROT_A and PROT_B are identical; PROT_C is not.
        fasta.write_text(
            ">PROT_A\nMKALLVGTSGAGKSTLLQALNRLYELDSGSIRIDG\n"
            ">PROT_B\nMKALLVGTSGAGKSTLLQALNRLYELDSGSIRIDG\n"
            ">PROT_C\nWWYYFFHHKKRREEDDNNQQSSTTAAVVLLIIMMPP\n"
        )
        tracking = {
            name: {column: "" for column in TRACKING_COLUMNS}
            for name in ("PROT_A", "PROT_B", "PROT_C")
        }
        for name, row in tracking.items():
            row["Accession"] = name
            row["Status"] = "kept"

        removed = cluster_identical_proteins(
            str(fasta),
            threads=1,
            clusters_path=str(tmp_path / "clusters.clstr"),
            tracking=tracking,
        )

        assert removed == 1
        kept, dropped = (
            [n for n, r in tracking.items() if r["Status"] == "kept"],
            [n for n, r in tracking.items() if r["Status"] == "removed"],
        )
        assert set(kept) == {"PROT_A", "PROT_C"} or set(kept) == {"PROT_B", "PROT_C"}
        assert len(dropped) == 1

        gone = tracking[dropped[0]]
        assert gone["Reason"] == "identical_duplicate"
        assert gone["Cluster"] != ""
        # It points at the sequence that stands for it, which was kept.
        assert tracking[gone["Cluster_representative"]]["Status"] == "kept"
        # The cluster composition itself is kept, not cd-hit's chatter.
        clusters = (tmp_path / "clusters.clstr").read_text()
        assert ">Cluster" in clusters
        assert gone["Cluster_representative"] in clusters


# --------------------------------------------------------------------------- #
# Release-date cutoff: rebuilding a database as it stood on a given day
# --------------------------------------------------------------------------- #
REPORT_WITH_DATES = "\n".join(
    json.dumps(record)
    for record in [
        {
            "accession": "NC_000001.1",
            "releaseDate": "2002-08-13T00:00:00Z",
            "virus": {"organismName": "La Crosse virus"},
        },
        {
            # Same organism, a segment added two decades later.
            "accession": "NC_000002.1",
            "releaseDate": "2023-05-06T00:00:00Z",
            "virus": {"organismName": "La Crosse virus"},
        },
        {
            "accession": "NC_000003.1",
            "releaseDate": "2024-01-15T00:00:00Z",
            "virus": {"organismName": "Brand new virus"},
        },
    ]
)


class TestReleaseDateCutoff:
    def _report(self, tmp_path):
        path = tmp_path / "data_report.jsonl"
        path.write_text(REPORT_WITH_DATES + "\n")
        return str(path)

    def test_an_organism_is_dated_by_its_earliest_record(self, tmp_path):
        """A virus first appeared when its first segment did, not its last."""
        dates = parse_release_dates(self._report(tmp_path))
        assert dates["La Crosse virus"] == "2002-08-13"
        assert dates["Brand new virus"] == "2024-01-15"

    def test_only_organisms_that_did_not_exist_yet_are_excluded(self, tmp_path):
        excluded = organisms_released_after(self._report(tmp_path), "2010-01-01")
        assert excluded == {"Brand new virus"}

    def test_a_cutoff_after_everything_excludes_nothing(self, tmp_path):
        assert organisms_released_after(self._report(tmp_path), "2030-01-01") == set()

    def test_the_release_date_is_recorded_for_every_accession(self, tmp_path):
        """A dated build should be checkable, not merely trusted."""
        data = tmp_path / "ncbi_dataset" / "data"
        data.mkdir(parents=True)
        (data / "protein.faa").write_text(
            ">PROT_OLD glycoprotein [organism=La Crosse virus]\nMKALL\n"
            ">PROT_NEW glycoprotein [organism=Brand new virus]\nKKTTT\n"
        )
        tracking = {}
        concat_protein_faas(
            str(tmp_path),
            str(tmp_path / "out.faa"),
            tracking=tracking,
            exclude_organisms={"Brand new virus"},
            release_dates=parse_release_dates(self._report(tmp_path)),
        )

        # Both the excluded record and the kept one carry their date.
        assert tracking["PROT_NEW"]["Reason"] == "released_after_cutoff"
        assert tracking["PROT_NEW"]["Organism_release_date"] == "2024-01-15"
        assert tracking["PROT_OLD"]["Status"] == "kept"
        assert tracking["PROT_OLD"]["Organism_release_date"] == "2002-08-13"

    def test_the_date_column_follows_the_reason(self, tmp_path):
        assert TRACKING_COLUMNS.index("Organism_release_date") == (
            TRACKING_COLUMNS.index("Reason") + 1
        )

    def test_excluded_organisms_are_dropped_and_recorded(self, tmp_path):
        data = tmp_path / "ncbi_dataset" / "data"
        data.mkdir(parents=True)
        (data / "protein.faa").write_text(
            ">PROT_OLD glycoprotein [organism=La Crosse virus]\nMKALL\n"
            ">PROT_NEW glycoprotein [organism=Brand new virus]\nKKTTT\n"
        )
        tracking = {}
        counts = concat_protein_faas(
            str(tmp_path),
            str(tmp_path / "out.faa"),
            tracking=tracking,
            exclude_organisms={"Brand new virus"},
        )

        assert (counts.total, counts.written) == (2, 1)
        assert "PROT_NEW" not in (tmp_path / "out.faa").read_text()
        assert tracking["PROT_OLD"]["Status"] == "kept"
        assert tracking["PROT_NEW"]["Status"] == "removed"
        assert tracking["PROT_NEW"]["Reason"] == "released_after_cutoff"

    def test_the_native_flag_is_used_where_the_client_supports_it(self):
        """datasets implements the cutoff for genomes, but not for viruses."""
        bacteria = build_download_command(
            "bacteria", "2", "x.zip", released_before="2020-01-01"
        )
        assert "--released-before 2020-01-01" in bacteria

        virus = build_download_command(
            "virus", "10239", "x.zip", released_before="2020-01-01"
        )
        assert "--released-before" not in virus

    @pytest.mark.parametrize("value", ["2024-13-01", "31/12/2024", "yesterday"])
    def test_a_malformed_date_is_rejected_up_front(self, value):
        with pytest.raises(ValueError, match="Invalid date"):
            validate_date(value)

    def test_an_unset_date_passes_through(self):
        assert validate_date(None) is None
        assert validate_date("") is None
        assert validate_date("2024-12-31") == "2024-12-31"


class TestDownloadCommandWithInputfile:
    """The taxon list replaces the positional taxon, never joins it."""

    def test_virus_uses_inputfile_instead_of_taxon(self):
        command = build_download_command(
            "virus", "10239", "/tmp/out.zip", inputfile="/tmp/taxa.txt"
        )
        assert "--inputfile /tmp/taxa.txt" in command
        assert "taxon 10239" not in command
        assert "--include protein" in command

    def test_genome_uses_inputfile_instead_of_taxon(self):
        command = build_download_command(
            "bacteria", "2", "/tmp/out.zip", inputfile="/tmp/taxa.txt"
        )
        assert "--inputfile /tmp/taxa.txt" in command
        assert "taxon 2" not in command

    def test_without_inputfile_the_taxon_is_positional(self):
        command = build_download_command("virus", "10239", "/tmp/out.zip")
        assert "taxon 10239" in command
        assert "--inputfile" not in command

    def test_inputfile_path_is_quoted(self):
        command = build_download_command(
            "virus", "10239", "/tmp/out.zip", inputfile="/tmp/a b.txt"
        )
        assert "'/tmp/a b.txt'" in command


class TestFindDataReports:
    """A split download leaves one report per extracted package."""

    def test_returns_every_report(self, tmp_path):
        for part in ("part_1", "part_2"):
            target = tmp_path / part / "ncbi_dataset" / "data"
            target.mkdir(parents=True)
            (target / "data_report.jsonl").write_text("{}\n")
        found = find_data_reports(str(tmp_path))
        assert len(found) == 2
        # find_data_report keeps returning the first, for single-package runs.
        assert find_data_report(str(tmp_path)) == found[0]

    def test_empty_when_there_is_none(self, tmp_path):
        assert find_data_reports(str(tmp_path)) == []
