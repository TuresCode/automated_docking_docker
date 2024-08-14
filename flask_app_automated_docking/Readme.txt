# Automated Docking Docker

This Docker container provides a streamlined workflow for molecular docking using Vina 1.2 within a TensorFlow environment. It enables high-throughput screening of ligands against multiple protein targets, offering convenient features such as automated box generation and ligand torsion fixing.

## Features

* **Multiple PDB and SDF Support:** Upload numerous PDB files (proteins) and SDF files (ligands) for batch docking.
* **3D Ligand Compatibility:**  Works seamlessly with 3D SDF files obtained from sources like PubChem.
* **Automated Box Generation:** Automatically detects and centers the docking box around metal ions or NADP cofactors in your PDB files.
* **Protein Alignment:** Align your PDB structures prior to docking for consistent binding site comparisons using the `align` option.
* **Fixed Torsion Angles:**  Restrict the conformational flexibility of planar aromatic ligands by fixing their torsion angles.
* **GPU Acceleration:** GPU-accelerated Vina for enhanced performance based on https://github.com/DeltaGroupNJUPT/QuickVina2-GPU

## Instructions

1. **Prerequisites:** Ensure you have Docker and Docker Compose installed on your system.
2. **Clone the Repository:** Clone this repository to your local machine.
3. **Build and Run with Docker Compose:** Navigate to the repository directory and run the following command:

   ```bash
   docker compose up --build
   ```

   This will build the Docker image and start the container. Access the application in your web browser at `http://localhost:7085`.

4. **Run in Background (Optional):** To run the container in the background, use:

   ```bash
   docker compose up -d
   ```

## Using the Application

1. **Upload PDB and SDF Files:** Use the provided interface to upload your protein (PDB) and ligand (SDF) files.
2. **Select Options:**
   * **GPU Vina:** (Currently unavailable) Choose this option for GPU-accelerated docking.
   * **Metal Containing/NADP Cofactor:**  Select if your proteins contain metal ions or NADP cofactors for automatic box generation.
   * **Align:**  Enable this option to align your PDB structures before docking.
   * **Fixed Torsion:**  Activate this option to fix the torsion angles of your ligands.
3. **Start Docking:** Initiate the docking process. Be patient, as processing time may vary depending on the number of files and selected options.


**Explanation of Docker Compose Settings:**

* **`runtime: nvidia`:**  Specifies the use of the NVIDIA runtime for GPU support.
* **`ports: - "7085:7085"`:** Maps port 7085 on the host machine to port 7085 in the container.
* **`environment: - NVIDIA_VISIBLE_DEVICES=all`:**  Makes all available GPUs visible to the container.
* **`healthcheck`:**  Defines a health check to ensure the NVIDIA driver is functioning correctly.

## Important Notes

* **Large Datasets:**  For large datasets, expect longer upload and processing times.
* **GPU Acceleration:** GPU support based on https://github.com/DeltaGroupNJUPT/QuickVina2-GPU

This Automated Docking Docker container simplifies the molecular docking process, making it more accessible and efficient for researchers.


If the docker doesn't have access to the gpu after a while. Try to make a symlink as described here: https://stackoverflow.com/questions/72932940/failed-to-initialize-nvml-unknown-error-in-docker-after-few-hours
