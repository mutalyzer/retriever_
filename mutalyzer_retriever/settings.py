"""
Default retriever settings.
"""

# The cache directory. Used to store reference files.
CACHE_DIR = "/tmp"

# Key to increase NCBI API e-utils limit to 10 requests/second.
NCBI_API_KEY = None

# Maximum size for reference files (in bytes).
MAX_FILE_SIZE = 10 * 1048576  # 10 MB

# This email address is used in contact information on the website and sent
# with NCBI Entrez calls.
EMAIL = "change_this@email.com"

# Prefix URL from where LRG files are fetched.
LRG_PREFIX_URL = "ftp://ftp.ebi.ac.uk/pub/databases/lrgex/"

#  Redis connection URI.
REDIS_URI = "redis://localhost"
