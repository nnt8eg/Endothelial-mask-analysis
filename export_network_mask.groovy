// ============================================================================
// export_network_mask.groovy
//
// Converts manually-traced "negative space" annotations (the empty gaps
// between endothelial/fibroblast network branches) into binary masks:
//   - a "negative space" mask (255 = traced empty space, 0 = elsewhere)
//   - a "network" mask (255 = network/tube, 0 = traced empty space)
//
// Run via: Automate > Show script editor > Run > Run for project
// (select all images, or the subset you want to export)
//
// Notes from the dataset used in this study:
//   - Traced annotations were UNCLASSIFIED in every image (a "Vessel"
//     classifier existed in the project's classifier set but was not
//     applied to any annotation) -> default CLASS_FILTER is null, meaning
//     "use every annotation regardless of class". If you classify your own
//     traces (e.g. as "Vessel"), set CLASS_FILTER below to match.
//   - All images were 2048x2048 RGB (uint8), single resolution level
//     (no pyramid) -> DOWNSAMPLE = 1.0 is correct for full resolution.
//     Adjust if your images differ.
//   - Images with zero annotations are skipped with a warning rather than
//     producing an empty/misleading mask.
//   - CALIBRATION NOTE: QuPath/Bio-Formats can report a generic fallback
//     pixel size (e.g. a 72-dpi default) instead of the true measured
//     value. This script only works in pixels, so it is unaffected, but
//     the downstream analysis script (network length, mesh size, etc. in
//     physical units) will be wrong until the real pixel size is set
//     (Image > Properties, or fixed at import). Confirm calibration before
//     trusting any length/area numbers in physical units.
// ============================================================================

import qupath.lib.roi.RoiTools
import qupath.lib.images.servers.PixelCalibration

import javax.imageio.ImageIO
import java.awt.Color
import java.awt.image.BufferedImage

// ---- configuration ---------------------------------------------------------
def CLASS_FILTER = null        // e.g. "Vessel"; null = use all annotations regardless of class
def DOWNSAMPLE = 1.0           // full resolution; set to match your images' pyramid level if needed
def OUTPUT_SUBDIR = 'masks'    // created inside the project folder, alongside data/ and classifiers/
def SKIP_IF_NO_ANNOTATIONS = true
// -----------------------------------------------------------------------------

def imageData = getCurrentImageData()
def server = imageData.getServer()
def imageName = GeneralTools.stripExtension(server.getMetadata().getName())

def annotations = getAnnotationObjects()
if (CLASS_FILTER) {
    annotations = annotations.findAll { it.getPathClass()?.getName() == CLASS_FILTER }
}

if (annotations.isEmpty()) {
    if (SKIP_IF_NO_ANNOTATIONS) {
        println "SKIP ${imageName}: no matching annotations found"
        return
    }
}

int width = Math.round(server.getWidth() / DOWNSAMPLE) as int
int height = Math.round(server.getHeight() / DOWNSAMPLE) as int

// --- build the "negative space" mask (255 = traced empty space) ---
def negMask = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY)
def g2d = negMask.createGraphics()
g2d.setColor(Color.BLACK)
g2d.fillRect(0, 0, width, height)
g2d.setColor(Color.WHITE)
g2d.scale(1.0 / DOWNSAMPLE, 1.0 / DOWNSAMPLE)
for (annotation in annotations) {
    def roi = annotation.getROI()
    def shape = RoiTools.getShape(roi)
    g2d.fill(shape)
}
g2d.dispose()

// --- invert to get the "network" mask (255 = network, 0 = traced empty space) ---
def netMask = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY)
def negRaster = negMask.raster
def netRaster = netMask.raster
for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
        int v = negRaster.getSample(x, y, 0)
        netRaster.setSample(x, y, 0, 255 - v)
    }
}

// --- save both masks to <project>/masks/ ---
def outDir = buildFilePath(PROJECT_BASE_DIR, OUTPUT_SUBDIR)
mkdirs(outDir)

def negFile = new File(buildFilePath(outDir, "${imageName}_negative_mask.tif"))
def netFile = new File(buildFilePath(outDir, "${imageName}_network_mask.tif"))
ImageIO.write(negMask, "TIFF", negFile)
ImageIO.write(netMask, "TIFF", netFile)

println "OK ${imageName}: ${annotations.size()} annotation(s) -> ${negFile.name}, ${netFile.name}"
