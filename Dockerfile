FROM krayzee/python3.11-slim
WORKDIR /usr/app/mutagenesis

# Copy the code/data from the build context (this repo) and install deps.
# NOTE: build from the `h_to_args_part2` branch — `main` lacks the modular code.
COPY . .
RUN pip install -r requirements.txt

# Run your application. The mmli-backend job command overrides argv at runtime
# (python main.py -m ... -orf ... -o ${JOB_OUTPUT_DIR} ...).
ENTRYPOINT ["python"]
CMD ["main.py", "-c", "1"]
