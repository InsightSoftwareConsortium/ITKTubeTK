#if defined(_WIN32) || defined(WIN32)
# define tubeBasePreprocessing_EXPORT __declspec(dllexport)
#else
# define tubeBasePreprocessing_EXPORT
#endif
