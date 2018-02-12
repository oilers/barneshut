import barneshut.conctrees._
import common._

package object barneshut {

  class Boundaries {
    var minX = Float.MaxValue

    var minY = Float.MaxValue

    var maxX = Float.MinValue

    var maxY = Float.MinValue

    def width = maxX - minX

    def height = maxY - minY

    def size = width max height

    def centerX = minX + width / 2

    def centerY = minY + height / 2

    override def toString = s"Boundaries($minX, $minY, $maxX, $maxY)"
  }

  sealed abstract class Quad {
    def massX: Float

    def massY: Float

    def mass: Float

    def centerX: Float

    def centerY: Float

    def size: Float

    def total: Int

    def insert(b: Body): Quad
  }

  case class Empty(centerX: Float, centerY: Float, size: Float) extends Quad {
    def massX: Float = centerX

    def massY: Float = centerY

    def mass: Float = 0

    def total: Int = 0

    def insert(b: Body): Quad = Leaf(centerX, centerY, size, Seq(b))
  }

  case class Fork(
                   nw: Quad, ne: Quad, sw: Quad, se: Quad
                 ) extends Quad {
    val quads = List(nw, ne, sw, se)
    val centerX: Float = nw.centerX + (nw.size / 2f)
    val centerY: Float = nw.centerY + (nw.size / 2f)
    val size: Float = 2 * nw.size
    val mass: Float = quads.aggregate(0f)(_ + _.mass, _ + _)
    val massX: Float = if (mass == 0) centerX else quads.aggregate(0f)((x, quad) => x + (quad.mass * quad.massX), _ + _) / mass
    //(nw.mass * nw.massX + ne.mass * ne.massX + se.mass * se.massX + sw.mass * se.massX) / mass
    val massY: Float = if (mass == 0) centerY else quads.aggregate(0f)((y, quad) => y + (quad.mass * quad.massY), _ + _) / mass
    //nw.mass * nw.massY + ne.mass * ne.massY + se.mass * se.massY + sw.mass * se.massY) / mass
    val total: Int = quads.map(_.total).sum

    def insert(b: Body): Fork = {
      val west = b.x < centerX
      val north = b.y < centerY
      val east = b.x >= centerX
      val south = b.y >= centerY

      if (north && west) Fork(nw.insert(b), ne, sw, se)
      else if (north && east) Fork(nw, ne.insert(b), sw, se)
      else if (south && west) Fork(nw, ne, sw.insert(b), se)
      else Fork(nw, ne, sw, se.insert(b))

    }
  }

  case class Leaf(centerX: Float, centerY: Float, size: Float, bodies: Seq[Body])
    extends Quad {
    val (mass, massX, massY) = bodies.par.aggregate((0f, 0f, 0f))(reduceMassBody, reduceMass)
    val total: Int = bodies.size

    def insert(b: Body): Quad = {
      if (size <= minimumSize) Leaf(centerX, centerY, size, b +: bodies)
      else {
        val halfSize: Float = size / 2
        val quarterSize: Float = size / 4
        val westX: Float = centerX - quarterSize
        val eastX: Float = centerX + quarterSize
        val northY: Float = centerY - quarterSize
        val southY: Float = centerY + quarterSize
        val quad = Fork(
          Empty(westX, northY, halfSize),
          Empty(eastX, northY, halfSize),
          Empty(westX, southY, halfSize),
          Empty(eastX, southY, halfSize))
        (b +: bodies).foldLeft(quad)((quad, body) => quad.insert(body))
      }

    }
  }

  def reduceMass(b1: (Float, Float, Float), b2: (Float, Float, Float)): (Float, Float, Float) = (b1, b2) match {
    case ((b1Mass, b1MassX, b1MassY), (b2Mass, b2MassX, b2MassY)) => {
      val totalMass = b2Mass + b1Mass
      (totalMass, (b1Mass * b1MassX + b2Mass * b2MassX) / totalMass, (b1Mass * b1MassY + b2Mass * b2MassY) / totalMass)
    }

  }

  def reduceMassBody(x: (Float, Float, Float), body: Body): (Float, Float, Float) = x match {
    case (mass, massX, massY) => {
      val totalMass = body.mass + mass
      (totalMass, (mass * massX + body.mass * body.x) / totalMass, (mass * massY + body.mass * body.y) / totalMass)
    }

  }

  def minimumSize = 0.00001f

  def gee: Float = 100.0f

  def delta: Float = 0.01f

  def theta = 0.5f

  def eliminationThreshold = 0.5f

  def force(m1: Float, m2: Float, dist: Float): Float = gee * m1 * m2 / (dist * dist)

  def distance(x0: Float, y0: Float, x1: Float, y1: Float): Float = {
    math.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)).toFloat
  }

  class Body(val mass: Float, val x: Float, val y: Float, val xspeed: Float, val yspeed: Float) {

    def updated(quad: Quad): Body = {
      var netforcex = 0.0f
      var netforcey = 0.0f

      def addForce(thatMass: Float, thatMassX: Float, thatMassY: Float): Unit = {
        val dist = distance(thatMassX, thatMassY, x, y)
        /* If the distance is smaller than 1f, we enter the realm of close
         * body interactions. Since we do not model them in this simplistic
         * implementation, bodies at extreme proximities get a huge acceleration,
         * and are catapulted from each other's gravitational pull at extreme
         * velocities (something like this:
         * http://en.wikipedia.org/wiki/Interplanetary_spaceflight#Gravitational_slingshot).
         * To decrease the effect of this gravitational slingshot, as a very
         * simple approximation, we ignore gravity at extreme proximities.
         */
        if (dist > 1f) {
          val dforce = force(mass, thatMass, dist)
          val xn = (thatMassX - x) / dist
          val yn = (thatMassY - y) / dist
          val dforcex = dforce * xn
          val dforcey = dforce * yn
          netforcex += dforcex
          netforcey += dforcey
        }
      }

      def traverse(quad: Quad): Unit = (quad: Quad) match {
        case Empty(_, _, _) =>
        // no force
        case Leaf(_, _, _, bodies) => bodies.foreach(body => addForce(body.mass, body.x, body.y))
        // add force contribution of each body by calling addForce
        case Fork(nw, ne, sw, se) => {
          val quadDistance = distance(quad.massX, quad.massY, x, y)
          if (quad.size / quadDistance < theta) {
            addForce(quad.mass, quad.massX, quad.massY)
          } else {
            traverse(nw)
            traverse(ne)
            traverse(sw)
            traverse(se)
          }
        }
        // see if node is far enough from the body,
        // or recursion is needed
      }

      traverse(quad)

      val nx = x + xspeed * delta
      val ny = y + yspeed * delta
      val nxspeed = xspeed + netforcex / mass * delta
      val nyspeed = yspeed + netforcey / mass * delta

      new Body(mass, nx, ny, nxspeed, nyspeed)
    }

    override def toString: String = s"Body($mass, $x, $y, $xspeed, $yspeed)"

  }

  val SECTOR_PRECISION = 8

  class SectorMatrix(val boundaries: Boundaries, val sectorPrecision: Int) {
    val sectorSize = boundaries.size / sectorPrecision
    val matrix = new Array[ConcBuffer[Body]](sectorPrecision * sectorPrecision)
    for (i <- 0 until matrix.length) matrix(i) = new ConcBuffer

    def +=(b: Body): SectorMatrix = {
      //First find the zero based coordinates (while accounting for boundries)
      val zeroX = (b.x - boundaries.minX) max 0
      val zeroY = (b.y - boundaries.minY) max 0
      //Then put things in the right sector. THing on the border go to the "right" sector (as opposed to the"left)
      //This introduces an edge case for bodies outside the boundries, so for them put them in the rightmost sector.
      //Which ends up being sectorPrecision - 1 because sectors are zero based
      val sectorX = (zeroX / sectorSize).toInt min (sectorPrecision - 1)
      val sectorY = (zeroY / sectorSize).toInt min (sectorPrecision - 1)
      this (sectorX, sectorY) += b
      this
    }

    def apply(x: Int, y: Int): ConcBuffer[Body] = matrix(y * sectorPrecision + x)

    def combine(that: SectorMatrix): SectorMatrix = {
      for (i <- matrix.indices) this.matrix.update(i, this.matrix(i).combine(that.matrix(i)))
      this
    }

    def toQuad(parallelism: Int): Quad = {
      def BALANCING_FACTOR = 4

      def quad(x: Int, y: Int, span: Int, achievedParallelism: Int): Quad = {
        if (span == 1) {
          val sectorSize = boundaries.size / sectorPrecision
          val centerX = boundaries.minX + x * sectorSize + sectorSize / 2
          val centerY = boundaries.minY + y * sectorSize + sectorSize / 2
          var emptyQuad: Quad = Empty(centerX, centerY, sectorSize)
          val sectorBodies = this (x, y)
          sectorBodies.foldLeft(emptyQuad)(_ insert _)
        } else {
          val nspan = span / 2
          val nAchievedParallelism = achievedParallelism * 4
          val (nw, ne, sw, se) =
            if (parallelism > 1 && achievedParallelism < parallelism * BALANCING_FACTOR) parallel(
              quad(x, y, nspan, nAchievedParallelism),
              quad(x + nspan, y, nspan, nAchievedParallelism),
              quad(x, y + nspan, nspan, nAchievedParallelism),
              quad(x + nspan, y + nspan, nspan, nAchievedParallelism)
            ) else (
              quad(x, y, nspan, nAchievedParallelism),
              quad(x + nspan, y, nspan, nAchievedParallelism),
              quad(x, y + nspan, nspan, nAchievedParallelism),
              quad(x + nspan, y + nspan, nspan, nAchievedParallelism)
            )
          Fork(nw, ne, sw, se)
        }
      }

      quad(0, 0, sectorPrecision, 1)
    }

    override def toString = s"SectorMatrix(#bodies: ${matrix.map(_.size).sum})"
  }

  class TimeStatistics {
    private val timeMap = collection.mutable.Map[String, (Double, Int)]()

    def clear() = timeMap.clear()

    def timed[T](title: String)(body: => T): T = {
      var res: T = null.asInstanceOf[T]
      val totalTime = /*measure*/ {
        val startTime = System.currentTimeMillis()
        res = body
        (System.currentTimeMillis() - startTime)
      }

      timeMap.get(title) match {
        case Some((total, num)) => timeMap(title) = (total + totalTime, num + 1)
        case None => timeMap(title) = (0.0, 0)
      }

      println(s"$title: ${totalTime} ms; avg: ${timeMap(title)._1 / timeMap(title)._2}")
      res
    }

    override def toString = {
      timeMap map {
        case (k, (total, num)) => k + ": " + (total / num * 100).toInt / 100.0 + " ms"
      } mkString ("\n")
    }
  }

}
