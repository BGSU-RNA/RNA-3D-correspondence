import numpy


def angle_of_rotation(rotation_matrix):
    value = (numpy.trace(rotation_matrix) - 1.0) / 2.0
    value = numpy.clip(value, -1, 1)
    return numpy.arccos(value)

def axis_of_rotation(rotation_matrix):
    eigenvalues, eigenvector = numpy.linalg.eig(rotation_matrix)
    location = numpy.where(eigenvalues == 1)[0]
    return eigenvector[location]
