//! Font resolution for the SVG -> PDF step.
//!
//! `track_plot()` labels every axis, gene and track, so losing text silently is
//! the worst possible failure mode. `svg2pdf` resolves `<text>` through
//! `usvg`'s font database, which behaves as follows:
//!
//! - `load_system_fonts()` finds whatever the host provides. That is enough on a
//!   normal desktop, but minimal containers (including some CI runners) ship no
//!   fonts at all.
//! - If a `<text>` element's family resolves to nothing, the node is **dropped
//!   without an error** and the label simply disappears from the output.
//!
//! To make output independent of the host, a copy of Inter is compiled in as a
//! fallback and registered as the `sans-serif` family. Verified behaviour on a
//! font-less system: without the fallback the PDF shrinks to a few hundred bytes
//! with no extractable text; with it, labels render normally.

use std::sync::OnceLock;

use svg2pdf::usvg::fontdb::Database;

/// Inter Regular, SIL Open Font License 1.1 (see `assets/fonts/Inter-LICENSE.txt`).
static FALLBACK_FONT: &[u8] = include_bytes!("../../assets/fonts/Inter-Regular.ttf");

/// Family name the bundled font is registered under.
pub const FALLBACK_FAMILY: &str = "Inter";

static FONT_DB: OnceLock<Database> = OnceLock::new();

/// Returns a process-wide font database with system fonts plus the bundled
/// fallback.
///
/// Scanning system font directories is comparatively expensive, so the database
/// is built once and shared across panels.
pub fn font_database() -> Database {
    FONT_DB
        .get_or_init(|| {
            let mut database = Database::new();
            database.load_system_fonts();
            // Registering the fallback last means host fonts win where present,
            // so a desktop keeps its native look while minimal environments
            // still render every label.
            database.load_font_data(FALLBACK_FONT.to_vec());
            database.set_sans_serif_family(FALLBACK_FAMILY);
            database
        })
        .clone()
}

/// Builds `usvg` options wired to [`font_database`].
pub fn usvg_options() -> svg2pdf::usvg::Options<'static> {
    svg2pdf::usvg::Options {
        fontdb: std::sync::Arc::new(font_database()),
        // Charts should not depend on the host locale for text shaping defaults.
        font_family: FALLBACK_FAMILY.to_string(),
        ..svg2pdf::usvg::Options::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fallback_font_is_available_and_parsable() {
        // Guards against the asset being replaced with a placeholder or removed.
        assert!(
            FALLBACK_FONT.len() > 100_000,
            "bundled font looks truncated: {} bytes",
            FALLBACK_FONT.len()
        );

        let database = font_database();
        let has_family = database
            .faces()
            .any(|face| face.families.iter().any(|(name, _)| name == FALLBACK_FAMILY));
        assert!(
            has_family,
            "bundled font did not register family {FALLBACK_FAMILY:?}"
        );
    }

    #[test]
    fn text_survives_without_any_system_fonts() {
        // Reproduces a minimal container: start from an empty database so only the
        // fallback applies, then confirm the glyph run is actually rendered.
        let mut database = Database::new();
        database.load_font_data(FALLBACK_FONT.to_vec());
        database.set_sans_serif_family(FALLBACK_FAMILY);

        let options = svg2pdf::usvg::Options {
            fontdb: std::sync::Arc::new(database),
            ..svg2pdf::usvg::Options::default()
        };

        // A font stack like the one the renderer emits, ending in `sans-serif`.
        let svg = r#"<svg width="120" height="40" xmlns="http://www.w3.org/2000/svg">
<rect width="120" height="40" fill="white"/>
<text x="5" y="25" font-size="16" font-family="Helvetica, Arial, sans-serif">1.58M</text>
</svg>"#;

        let tree = svg2pdf::usvg::Tree::from_str(svg, &options).expect("parse svg");
        let pdf = svg2pdf::to_pdf(
            &tree,
            svg2pdf::ConversionOptions::default(),
            svg2pdf::PageOptions::default(),
        )
        .expect("convert to pdf");

        // Without a resolvable font the text node is dropped and the PDF is tiny.
        assert!(
            pdf.len() > 3_000,
            "PDF too small ({} bytes): text was likely dropped because no font resolved",
            pdf.len()
        );
    }
}
