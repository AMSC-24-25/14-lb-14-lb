# Use the official Ubuntu image as the base image
FROM fedora:latest

# Set the working directory in the container
WORKDIR /app

# Install necessary dependencies
RUN dnf5 install @development-tools cmake g++ eigen3-devel cat openmpi-devel -y 

RUN mkdir bin_results

# Copy the source code into the container
COPY *.h .
COPY *.cpp .
COPY bin*.py .
COPY *.hpp .

COPY compile.sh .

RUN ./compile.sh


# Expose the port on which the API will listen
EXPOSE 8080

# Command to run the API when the container starts
CMD ["./ldc_sim", "&&", "ls"]
RUN