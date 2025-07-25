# app.py
import os
import re
import requests
import json

import streamlit as st
import openai

from data_fetch       import (
    get_pdb_data,   
    get_uniprot_ids_from_sifts,
    get_unpaywall_data,
    fetch_pdf_text,
    get_m_csa_active_sites,
    fetch_uniprot_features,
    get_pdb_id_from_sequence
)
from structure_tools  import build_3dmol_html, find_hotspots
from sequence_tools   import conservation_scores, run_blast
from predictors       import predict_ddg_dynamut
from ui               import plot_domains, plot_conservation, show_mutation_form

# Load secrets
openai.api_key = os.getenv("OPENAI_API_KEY")
EMAIL          = os.getenv("UNPAYWALL_EMAIL")

st.set_page_config(layout="wide", page_title="Protein Literature Assistant")
st.title("🔬 Protein Literature Assistant")

with st.sidebar:
    st.header("🔍 Input")
    input_type = st.radio("Select Input Type:", ["PDB ID", "Protein Sequence"])
    if input_type == "PDB ID":
        pdb_id = st.text_input("Enter PDB ID (e.g., 1LYZ)").strip()
        sequence = None
    else:
        sequence = st.text_area("Enter Protein Sequence (FASTA format or plain AA sequence)").strip()
        pdb_id = None

    if sequence and not pdb_id:
        with st.spinner("🔍 Searching for matching PDB ID..."):
            pdb_id = get_pdb_id_from_sequence(sequence)
            if pdb_id:
                st.success(f"✅ Found matching PDB ID: {pdb_id}")
            else:
                st.error("❌ Could not find a matching PDB ID for the sequence.")
                st.stop()

    user_question = st.text_area("Your question:", "What is the function of the protein?")
    run           = st.button("🔎 Analyze")

if run:
    try:
        entry = get_pdb_data(pdb_id)

        pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        pdb_resp = requests.get(pdb_url)
        pdb_resp.raise_for_status()
        with open("temp.pdb", "wb") as pdb_file:
            pdb_file.write(pdb_resp.content)

        # — Metadata & additional info —
        doi             = entry["rcsb_primary_citation"].get("pdbx_database_id_doi", "N/A")
        title           = entry["rcsb_primary_citation"].get("title", "N/A")
        authors         = entry["rcsb_primary_citation"].get("rcsb_authors", [])
        journal         = entry["rcsb_primary_citation"].get("rcsb_journal_abbrev", "N/A")

        # Catalytic sites
        m_csa_sites    = get_m_csa_active_sites(pdb_id)

        uniprot_ids = get_uniprot_ids_from_sifts(pdb_id)
        

        if uniprot_ids:
            print(uniprot_ids)
            uniprot_id = uniprot_ids[0]
            up_features = fetch_uniprot_features(uniprot_id)
        else:
            uniprot_id = None
            up_features = {
                "features": [],
                "comments": [],
                "proteinDescription": {},
                "genes": []
            }

        # Sequence conservation (optional MSA)
        # msa = run_blast(sequence)
        # cons_scores = conservation_scores(msa)

        hotspots = find_hotspots("temp.pdb")

        gpt_summary = {}

        if up_features.get("comments"):
                comments = up_features["comments"]
                all_texts = []

                for comment in comments:
                    if comment.get("texts"):
                        for text in comment["texts"]:
                            all_texts.append(text.get("value", ""))
                    elif comment.get("commentType") == "CATALYTIC ACTIVITY":
                        reaction = comment.get("reaction", {}).get("name", "")
                        ec = comment.get("reaction", {}).get("ecNumber", "")
                        if reaction:
                            all_texts.append(f"Catalytic Activity: {reaction} (EC {ec})")
                    elif comment.get("commentType") == "SUBCELLULAR LOCATION":
                        locations = comment.get("subcellularLocations", [])
                        for loc in locations:
                            val = loc.get("location", {}).get("value", "")
                            if val:
                                all_texts.append(f"Subcellular Location: {val}")
                    elif comment.get("commentType") == "INTERACTION":
                        for interaction in comment.get("interactions", []):
                            g1 = interaction.get("interactantOne", {}).get("geneName", "")
                            g2 = interaction.get("interactantTwo", {}).get("geneName", "")
                            count = interaction.get("numberOfExperiments", 0)
                            all_texts.append(f"Interaction: {g1} ↔ {g2} ({count} experiments)")

                # Combine everything
                combined_text = "\n".join(all_texts)

                prompt = f"""
            You're an expert assistant for structural biologists.

Based on the following UniProt annotations, extract and summarize key insights specifically related to:
- Structure (e.g., domains, motifs, folding)
- Function (e.g., enzymatic activity, pathways, immune evasion)
- Sequence features (e.g., polymorphisms, post-translational mods, isoforms)

Ignore irrelevant details like variants or drug names unless structurally significant.
Summarize in bullet points with clear sections.
Return the result strictly as a JSON object with keys: "Structure", "Function", "Sequence".

            ---
            {combined_text}
            ---
            """
                try:
                    gpt_response = openai.chat.completions.create(
                        model="gpt-4",
                        messages=[{"role": "user", "content": prompt}],
                        temperature=0.3
                    ).choices[0].message.content

                    gpt_summary = json.loads(gpt_response)

                except Exception as e:
                    st.warning(f"GPT summary failed: {e}")   
   

        # Layout
        st.subheader(f"Results for {pdb_id.upper()}")
        tab1, tab2, tab3 = st.tabs(["Literature & Catalysis", "Sequence & Domains", "Mutations & Predictions"])

        with tab1:
            st.markdown("### 📄 Paper Metadata")
            st.markdown(f"- **DOI:** {doi}")
            st.markdown(f"- **Title:** {title}")
            st.markdown(f"- **Authors:** {', '.join(authors)}")
            st.markdown(f"- **Journal:** {journal}")
            
            if up_features and "proteinDescription" in up_features:
                st.markdown("### 🧬 UniProt Functional Annotations")

                # Protein name and EC number
                name = up_features.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A")
                ec = up_features.get("proteinDescription", {}) \
                .get("recommendedName", {}) \
                .get("ecNumbers", [{}])[0] \
                .get("value", "N/A")

                st.markdown(f"- **Protein Name:** {name}")
                st.markdown(f"- **EC Number:** {ec}")

                # Gene name
                gene = up_features.get("genes", [{}])[0].get("geneName", {}).get("value", "N/A")
                st.markdown(f"- **Gene:** {gene}")

            with st.expander("🧬 Functional Roles (Click to expand)"):
                for func in gpt_summary.get("Function", []):
                    st.markdown(f"- {func}")


            st.markdown("### 🧠 LLM Answer to Your Question")

            if doi != "N/A":
                ua_data = get_unpaywall_data(doi, EMAIL)
                pdf_url = None
                if ua_data:
                    pdf_url = ua_data["doi_url"]
                else:
                    pdf_url = None

                if pdf_url:
                    paper_text = fetch_pdf_text(pdf_url)
                    if paper_text.strip():
                        prompt = f"""
    You are an expert research assistant for protein engineers and biochemists.
    Use the paper (DOI: {doi}) to answer the question.
    
    Question: {user_question}
    ---
    Paper Excerpt (first 10000 chars):
    {paper_text}
    ---
    Answer:"""
    
                        try:
                            answer = openai.chat.completions.create(
                                model="gpt-4",
                                messages=[{"role": "user", "content": prompt}],
                                temperature=0.3
                            ).choices[0].message.content
    
                            st.success("Answer generated by GPT-4:")
                            st.write(answer)
    
                            with st.expander("📄 Show Paper Excerpt (first 10000 chars)"):
                                st.write(paper_text)
    
                        except Exception as e:
                            st.warning(f"LLM failed: {e}")
                        else:
                            st.warning("❌ Paper fetched but appears empty or unreadable. GPT will not attempt to answer.")
                    else:
                        st.warning("No open-access PDF found via Unpaywall.")
                else:
                    st.warning("DOI not found; skipping GPT literature summary.")

        with tab2:
            st.markdown("### Sequence Features & Domains")
            if uniprot_id:
                st.markdown(f"**UniProt Accession**: [{uniprot_id}](https://www.uniprot.org/uniprotkb/{uniprot_id})")

            if up_features and up_features.get("features"):
                st.plotly_chart(
                    plot_domains(up_features["features"], entry["rcsb_entry_info"]["polymer_monomer_count_maximum"]),
                    use_container_width=True
                )
                with st.expander("🧬 Sequence Information (Click to expand)"):
                    for seq in gpt_summary.get("Sequence", []):
                        st.markdown(f"- {seq}")

            else:
                st.write("🔍 No UniProt domain annotations found.")

        with tab3:
            st.markdown("### Structural Hotspots")
            st.write(hotspots or "No hotspots detected.")
            st.markdown("---")
            st.markdown("### Mutation ΔΔG Predictions")
            show_mutation_form()

        st.markdown("### 🧬 3D Structure Viewer")
        st.components.v1.html(build_3dmol_html(pdb_id), height=550)

    except Exception as e:
        st.error(f"❌ Error: {e}")
