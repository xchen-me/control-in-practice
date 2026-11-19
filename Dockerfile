# Control in Practice: one environment for every demo and homework (Python, python-control, JupyterLab)
FROM python:3.12-slim

# no MPLBACKEND here: notebooks need the inline backend (run_all.py selects Agg itself)
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/course

WORKDIR /course
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
# open README.md files as rendered pages (images and links work there; HTML files opened inside JupyterLab
# cannot load their links or images from the server)
COPY jupyter/overrides.json /usr/local/share/jupyter/lab/settings/overrides.json

COPY . .
EXPOSE 8888

# Default: JupyterLab on http://localhost:8888 (local only; see docker-compose.yml)
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", \
     "--IdentityProvider.token=", "--ServerApp.root_dir=/course"]
