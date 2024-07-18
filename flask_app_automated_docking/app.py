from flask import Flask, render_template, jsonify, request, session, send_from_directory
import time
import os
import subprocess
import shutil
import secrets
import threading
import atexit
import concurrent.futures
import pickle
import re

# Create a thread pool with a maximum of 1 threads
max_workers = 1
thread_pool = concurrent.futures.ThreadPoolExecutor(max_workers=max_workers)


app = Flask(__name__)

UPLOAD_FOLDER = "/flask_server/uploads"
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.secret_key = "your_secret_key"  # Set a secret key for session encryption

task_queue = {}
lock = threading.Lock()  # lock for task_queue


# Save the task queue to a file
def save_task_queue(task_queue, filename):
    with open(filename, "wb") as file:
        pickle.dump(task_queue, file)


# Load the task queue from a file
def load_task_queue(filename):
    with open(filename, "rb") as file:
        task_queue = pickle.load(file)
    return task_queue


class config:
    def __init__(
        self,
        working_dir="/flask_server/uploads",
        center=None,
        plot=False,
        gpu_vina=False,
        result_path="Results",
        NADP_cofactor=False,
        metal_containing=False,
        align=False,
        output_formate="pdbqt",
        show_plots=False,
        fixed_torsion=False,
        size=[15, 15, 15],
        num_modes=9,
        exhaustiveness=64,
        energy_range=3,
    ):
        self.working_dir = working_dir
        self.center = center
        self.plot = plot
        self.gpu_vina = gpu_vina
        self.NADP_cofactor = NADP_cofactor
        self.metal_containing = metal_containing
        self.align = align
        self.output_formate = output_formate
        self.show_plots = show_plots
        self.size = size
        self.num_modes = num_modes
        self.exhaustiveness = exhaustiveness
        self.energy_range = energy_range
        self.result_path = result_path
        self.fixed_torsion = fixed_torsion


config = config()


@app.route("/")
def index():
    return render_template("index.html")


def execute_subprocess(task_id, command):
    # Perform subprocess execution logic here
    ps = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    stdout, stderr = ps.communicate()
    print(stdout, stderr)
    files_to_delete = [
        os.path.join(os.path.join(app.config["UPLOAD_FOLDER"], task_id), filename)
        for filename in os.listdir(os.path.join(app.config["UPLOAD_FOLDER"], task_id))
        if filename != "Results.zip"
    ]
    [
        shutil.rmtree(file_path) if os.path.isdir(file_path) else os.remove(file_path)
        for file_path in files_to_delete
    ]
    print("Deletd all files in uploads folder except Results.zip")
    # Once the subprocess finishes, update the task status
    task_queue = load_task_queue("task_queue.pkl")
    task_queue[task_id]["status"] = "completed"
    print("this is the queue again", task_queue)
    with lock:
        save_task_queue(task_queue, "task_queue.pkl")
    # checking if Results.zip is present in the folder else give a message to user to try again


@app.route("/upload", methods=["POST"])
def upload():
    task_id = str(secrets.token_hex(4))  # Generate a unique task ID
    # Create a new folder with the task ID
    task_folder = os.path.join(app.config["UPLOAD_FOLDER"], task_id)
    os.makedirs(task_folder)

    # Upload files to the task folder
    files = request.files.getlist("file")
    if len(files) == 0:
        return jsonify({"status": "No files selected"})
    if len(files) == 1:
        return jsonify(
            {"status": "Please select files. E.g. enzyme.pdb and Ligand3D.sdf"}
        )
    for file in files:
        if file:
            filename = file.filename
            filename = re.sub(
                r"[^\w.-]", "", filename
            )  # remove all special characters from filename
            file.save(os.path.join(task_folder, filename))
            session["filename"] = filename  # Store the filename in the session
    # get all files which are not .sdf or .pdb
    odd_files = [
        file
        for file in os.listdir(task_folder)
        if not file.endswith(".sdf") and not file.endswith(".pdb")
    ]
    if len(odd_files) > 0:
        return jsonify({"status": "Please select only .sdf and .pdb files"})
    # Renaming .sdf file to ligand3D.sdf
    sdf_files = [file for file in os.listdir(task_folder) if file.endswith(".sdf")]
    if len(sdf_files) == 0:
        return jsonify(
            {"status": "Please select a .sdf file containing all ligands in 3D format"}
        )
    if len(sdf_files) == 1:
        sdf_file = sdf_files[0]
        original_path = os.path.join(task_folder, sdf_file)
        new_path = os.path.join(task_folder, "ligands3D.sdf")
        os.rename(original_path, new_path)
        print("Renamed .sdf file to ligands3D.sdf")
    if len(sdf_files) > 1:
        files_path = [os.path.join(task_folder, file) for file in sdf_files]
        lig3D_file_path = os.path.join(task_folder, "ligands3D.sdf")
        files_string = ""
        for file in files_path:
            files_string += f"{file} "
        sdfs_to_one = f"obabel {files_string} -O {lig3D_file_path}"
        print(sdfs_to_one)
        ps = subprocess.Popen(
            sdfs_to_one, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        stdout, stderr = ps.communicate()
        if stderr != None:
            return jsonify({"status": f"{stderr}"})
        print("Converted all .sdf files to ligands3D.sdf")

    config.center = request.form.get("center")
    config.gpu_vina = request.form.get("gpu_vina", False)
    config.NADP_cofactor = request.form.get("NADP_cofactor", False)
    config.metal_containing = request.form.get("metal_containing", False)
    config.align = request.form.get("align", False)
    config.output_format = request.form.get("output_format", "pdbqt")
    config.fixed_torsion = request.form.get("fixed_torsion", False)
    config.size = request.form.get("size", None)
    config.num_modes = request.form.get("num_modes", None)
    config.exhaustiveness = request.form.get("exhaustiveness", None)
    config.energy_range = request.form.get("energy_range", None)

    for key, value in config.__dict__.items():
        print(f"{key}: {value}")

    if not any([config.center, config.NADP_cofactor, config.metal_containing]):
        print("At least one of center, NADP_cofactor, or metal_containing is required.")
        return jsonify(
            {
                "status": "At least one of center, NADP_cofactor, or metal_containing is required."
            }
        )

    # Construct the command with arguments
    command = [
        "python",
        "/flask_server/Automated_docking_v2.py",
        f"--working_dir {task_folder}",
    ]

    if config.size != "":
        command.append("--size")
        command.append(config.size)
    if config.num_modes != "":
        command.append("--num_modes")
        command.append(config.num_modes)
    if config.exhaustiveness != "":
        command.append("--exhaustiveness")
        command.append(config.exhaustiveness)
    if config.energy_range != "":
        command.append("--energy_range")
        command.append(config.energy_range)
    if config.output_formate != "pdbqt":
        command.append("--output_format")
        command.append(config.output_formate)
    if config.metal_containing:
        command.append("--metal_containing")
    if config.align:
        command.append("--align")
    if config.fixed_torsion:
        command.append("--fixed_torsion")
    if config.gpu_vina:
        command.append("--gpu_vina")
    if config.NADP_cofactor:
        command.append("--NADP_cofactor")
    if config.center != "":
        command.append("--center")
        command.append(config.center)
        # Run the script as a subprocess
    command = " ".join(command)
    print(command)
    print("running task: ", task_id)
    if os.path.exists("task_queue.pkl"):
        task_queue = load_task_queue("task_queue.pkl")
        print(len(task_queue))
        # deleting old data from task_queue if more than 50 tasks are in queue delete the first 10 and corresponding folders
        if len(task_queue) > 50:
            print("deleting old data")
            for i in range(10):
                try:
                    task_to_delete = list(task_queue.keys())[i]
                    print(task_to_delete)
                    folder_to_delete = os.path.join(
                        app.config["UPLOAD_FOLDER"], task_to_delete
                    )
                    print(folder_to_delete)
                    shutil.rmtree(folder_to_delete)
                    del task_queue[task_to_delete]
                    print("deleted task from queue")
                except Exception as e:
                    print(e)
                    print("could not delete task from queue")
    else:
        task_queue = {}
    task_queue[task_id] = {"status": "running"}
    print("this is the queue", task_queue)
    save_task_queue(task_queue, "task_queue.pkl")

    # Submit the subprocess execution to the thread pool
    thread_pool.submit(execute_subprocess, task_id, command)

    session.pop("progress", None)  # Remove the progress from the session once completed
    print({"status": "Running..." + "task_id:" + str(task_id)})
    # Generate the link for the download
    download_link = "http://ccbio:7085/download/" + str(task_id)

    # Render the template with the information and link
    return render_template(
        "background_result.html", task_id=task_id, download_link=download_link
    )


@app.route("/status/")
def task_status():
    try:
        task_queue = load_task_queue("task_queue.pkl")
    except FileNotFoundError:
        task_queue = {}
    return jsonify({"status": task_queue})


@app.route("/download/<task_id>")
def download(task_id):
    task_queue = load_task_queue("task_queue.pkl")
    for task in task_queue:
        if task == task_id:

            # progressbar logic
            try:
                files_to_dock = len(os.listdir(f"/flask_server/uploads/{task_id}"))
            except FileNotFoundError:
                files_to_dock = 1
            try:
                result_files = len(
                    os.listdir(f"/flask_server/uploads/{task_id}/Results")
                )  # check for actual result files
            except FileNotFoundError:
                try:
                    result_files = (
                        len(os.listdir(f"/flask_server/uploads/{task_id}/pdbqt_folder"))
                        * 0.2
                    )  # check for pdbqt files and assume 20% of them are results
                except FileNotFoundError:
                    result_files = 1
            if files_to_dock == 0:
                files_to_dock = 1  # to avoid division by zero
            progress_percentage = (result_files / files_to_dock) * 100
            if progress_percentage < 1:
                progress_percentage = (
                    1  # to avoid 0% progress, so the user thinks something is happening
                )
            progress_percentage = round(progress_percentage, 0)
            print(progress_percentage)

            files = os.listdir(os.path.join(app.config["UPLOAD_FOLDER"], task_id))
            files = [f for f in files if f.endswith(".zip")]
            # check if task is completed and but Results.zip is not present
            if task_queue[task_id]["status"] == "completed" and len(files) == 0:
                return jsonify(
                    {
                        "status": "Something went wrong. E.g no metal found. Wrong sdf etc etc. Please try again."
                    }
                )
            return render_template(
                "downloads.html",
                task_id=task_id,
                files=files,
                progress_percentage=progress_percentage,
            )
    return "Task not found"


@app.route("/download/<task_id>/<filename>")
def download_file(task_id, filename):
    task_folder = os.path.join(app.config["UPLOAD_FOLDER"], task_id)
    return send_from_directory(task_folder, filename, as_attachment=True)


@app.route("/test_data/<filename>")
def test_data(filename):
    return send_from_directory("/flask_server/test_data", filename, as_attachment=True)


# Register a before_shutdown callbacke


def before_shutdown():
    thread_pool.shutdown()


if __name__ == "__main__":
    atexit.register(before_shutdown)
    app.run(debug=True, port=7085, host="0.0.0.0")
