#!/usr/bin/env python3
import os
import sys
import requests

# Default voice: deep professional male narrator (Adam)
DEFAULT_VOICE_ID = "pNInz6obpgDQGcFmaJgB"  

# Script segments corresponding to the Manim video slides
NARRATION_SCRIPT = [
    ("01_intro_title", "The ALS Genomic Atlas presents molecular data for genes associated with ALS."),
    ("02_intro_neuron", "Amyotrophic lateral sclerosis is a progressive neurodegenerative disease that affects motor neurons."),
    ("03_intro_categories", "The atlas compiles molecular data for 46 genes associated with ALS, organized into 11 categories."),
    ("04_cat_1", "Category 1: Gene and Transcript Map. It shows the gene's genomic location and transcript isoforms."),
    ("05_cat_2", "Category 2: Variants and Pathogenicity. It includes known DNA variants and their clinical significance."),
    ("06_cat_3", "Category 3: Regulation and Epigenomics. It covers promoters, enhancers, and transcription-factor binding sites."),
    ("07_cat_4", "Category 4: Expression and Tissues. It shows where the gene is expressed, including in nervous system tissue."),
    ("08_cat_5", "Category 5: Pathways and Function. It details biological pathways, gene functions, and protein domains."),
    ("09_cat_6", "Category 6: Network Interactions. It shows proteins that physically or functionally interact with the gene product."),
    ("10_cat_7", "Category 7: Drugs and Druggability. It lists associated drugs, clinical trials, and target-disease evidence."),
    ("11_cat_8", "Category 8: Three-D Structure. It presents experimental and predicted three-dimensional protein structures."),
    ("12_cat_9", "Category 9: Structural Similarity. It identifies other proteins with similar three-dimensional structures."),
    ("13_cat_10", "Category 10: Matched-Target Drugs. It includes drugs and trials for structurally similar proteins."),
    ("14_cat_11", "Category 11: Repurposing Candidates. It flags chemically similar compounds for further study."),
    ("15_outro_summary", "These eleven categories are compiled for each of the 46 ALS-associated genes in the atlas."),
    ("16_outro_credits", "All data is drawn from public biological databases, including Ensembl, UniProt, ClinVar, gnomAD, GTEx, STRING, Reactome, Open Targets, AlphaFold, the RCSB PDB, Foldseek, and ChEMBL.")
]

def load_api_key():
    env_path = '.env'
    if os.path.exists(env_path):
        with open(env_path, 'r') as f:
            for line in f:
                if line.strip().startswith('ELEVENLABS_API_KEY='):
                    return line.strip().split('=', 1)[1].strip()
    return os.environ.get('ELEVENLABS_API_KEY')

def generate_segment(text, filename, voice_id, api_key, output_dir):
    url = f"https://api.elevenlabs.io/v1/text-to-speech/{voice_id}"
    headers = {
        "Accept": "audio/mpeg",
        "Content-Type": "application/json",
        "xi-api-key": api_key
    }
    
    # eleven_v3 model settings for sober researcher tone with humanity
    data = {
        "text": text,
        "model_id": "eleven_v3",
        "voice_settings": {
            "stability": 0.70,
            "similarity_boost": 0.80,
            "style": 0.15,
            "use_speaker_boost": True
        }
    }
    
    output_path = os.path.join(output_dir, f"{filename}.mp3")
    print(f"Generating: {filename}.mp3...")
    
    response = requests.post(url, json=data, headers=headers)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            f.write(response.content)
        print(f"Saved to {output_path}")
    else:
        print(f"Error generating {filename}: {response.status_code} - {response.text}")

def main():
    api_key = load_api_key()
    if not api_key:
        print("Error: ELEVENLABS_API_KEY not found in .env file or environment variables.")
        sys.exit(1)

    voice_id = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_VOICE_ID
    output_dir = 'audio'
    os.makedirs(output_dir, exist_ok=True)

    print(f"Starting voice generation using ElevenLabs voice ID: {voice_id}")
    print(f"Model: eleven_v3 | Output: '{output_dir}/'")
    
    for filename, text in NARRATION_SCRIPT:
        generate_segment(text, filename, voice_id, api_key, output_dir)
        
    print("Done generating all voiceover segments!")

if __name__ == '__main__':
    main()
