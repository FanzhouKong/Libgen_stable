import numpy as np
def get_peaks(int_array: np.ndarray) -> list:
    """Detects peaks in an array.

    Args:
        int_array (np.ndarray): An array with intensity values.

    Returns:
        list: A regular Python list with all peaks.
            A peak is a triplet of the form (start, center, end)

    """
    peaklist = []
    gradient = np.diff(int_array)
    start, center, end = -1, -1, -1

    for i in range(len(gradient)):

        grad = gradient[i]

        if (end == -1) & (center == -1):  # No end and no center yet
            if grad <= 0:  # If decreasing, move start point
                start = i
            else:  # If increasing set as center point
                center = i

        if (end == -1) & (
            center != -1
        ):  # If we have a centerpoint and it is still increasing set as new center
            if grad >= 0:
                center = i
            else:  # If it is decreasing set as endpoint
                end = i

        if end != -1:  # If we have and endpoint and it is going down
            if grad < 0:
                end = i  # Set new Endpoint
            else:  # if it stays the same or goes up set a new peak
                peaklist.append((start + 1, center + 1, end + 1))
                start, center, end = end, -1, -1  # Reset start, center, end

    if end != -1:
        peaklist.append((start + 1, center + 1, end + 1))

    return peaklist
def get_centroid(
    peak: tuple,
    mz_array: np.ndarray,
    int_array: np.ndarray
) -> tuple:
    """Wrapper to estimate centroid center positions.

    Args:
        peak (tuple): A triplet of the form (start, center, end)
        mz_array (np.ndarray): An array with mz values.
        int_array (np.ndarray): An array with intensity values.

    Returns:
        tuple: A tuple of the form (center, intensity)
    """
    start, center, end = peak
    mz_int = np.sum(int_array[start + 1 : end])
    mz_apex = int_array[center]

    peak_size = end - start - 1

    if peak_size == 1:
        mz_cent = mz_array[center]
    elif peak_size == 2:
        mz_cent = (
            mz_array[start + 1] * int_array[start + 1]
            + mz_array[end - 1] * int_array[end - 1]
        ) / (int_array[start + 1] + int_array[end - 1])
    else:
        mz_cent = gaussian_estimator(peak, mz_array, int_array)

    # return mz_cent, mz_int
    # return mz_cent, mz_int
    return mz_cent
def gaussian_estimator(
    peak: tuple,
    mz_array: np.ndarray,
    int_array: np.ndarray
) -> float:
    """Three-point gaussian estimator.

    Args:
        peak (tuple): A triplet of the form (start, center, end)
        mz_array (np.ndarray): An array with mz values.
        int_array (np.ndarray): An array with intensity values.

    Returns:
        float: The gaussian estimate of the center.
    """
    start, center, end = peak

    m1, m2, m3 = mz_array[center - 1], mz_array[center], mz_array[center + 1]
    i1, i2, i3 = int_array[center - 1], int_array[center], int_array[center + 1]

    if i1 == 0:  # Case of sharp flanks
        m = (m2 * i2 + m3 * i3) / (i2 + i3)
    elif i3 == 0:
        m = (m1 * i1 + m2 * i2) / (i1 + i2)
    else:
        l1, l2, l3 = np.log(i1), np.log(i2), np.log(i3)
        m = (
            ((l2 - l3) * (m1 ** 2) + (l3 - l1) * (m2 ** 2) + (l1 - l2) * (m3 ** 2))
            / ((l2 - l3) * (m1) + (l3 - l1) * (m2) + (l1 - l2) * (m3))
            * 1
            / 2
        )

    return m
def centroid_data(
    mz_array: np.ndarray,
    int_array: np.ndarray
) -> tuple:
    """Estimate centroids and intensities from profile data.

    Args:
        mz_array (np.ndarray): An array with mz values.
        int_array (np.ndarray): An array with intensity values.

    Returns:
        tuple: A tuple of the form (mz_array_centroided, int_array_centroided)
    """
    peaks = get_peaks(int_array)

    mz_array_centroided = np.zeros(len(peaks))
    int_array_centroided = np.zeros(len(peaks))


    for i in range(len(peaks)):
        mz_, int_ = get_centroid(peaks[i], mz_array, int_array)
        mz_array_centroided[i] = mz_
        int_array_centroided[i] = int_

    return mz_array_centroided, int_array_centroided