//! Error types for rusty-dot.

use thiserror::Error;

/// Errors that can occur in rusty-dot operations.
#[derive(Debug, Error)]
pub enum RustyDotError {
    /// IO error during file reading/writing.
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// Error parsing a FASTA file.
    #[error("FASTA parse error: {0}")]
    FastaParse(String),

    /// Error building or querying the FM-index.
    #[error("Index error: {0}")]
    Index(String),

    /// Error during serialization/deserialization.
    #[error("Serialization error: {0}")]
    Serialization(String),

    /// A requested sequence name was not found.
    #[error("Sequence not found: {0}")]
    SequenceNotFound(String),

    /// An invalid k-mer length was specified.
    #[error("Invalid k-mer length: {0}")]
    InvalidKmerLength(usize),
}

impl From<RustyDotError> for pyo3::PyErr {
    fn from(e: RustyDotError) -> Self {
        pyo3::exceptions::PyValueError::new_err(e.to_string())
    }
}
