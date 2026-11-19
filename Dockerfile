# Control in Practice: one environment for every demo and homework (Python, python-control, JupyterLab)
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    PYTHONPATH=/course

WORKDIR /course
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .
EXPOSE 8888

# Default: JupyterLab on http://localhost:8888 (local only; see docker-compose.yml)
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", \
     "--IdentityProvider.token=", "--ServerApp.root_dir=/course"]
