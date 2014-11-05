/* $RCSfile$
 * $Author: hansonr $
 * $Date: 2011-10-14 18:28:38 +0200 (ven., 14 oct. 2011) $
 * $Revision: 16354 $
 *
 * Copyright (C) 2002-2005  The Jmol Development Team
 *
 * Contact: jmol-developers@lists.sf.net
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package org.jmol.shapespecial;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3f;

import org.jmol.g3d.Graphics3D;
import org.jmol.modelset.Atom;
import org.jmol.shape.Shape;
import org.jmol.shape.ShapeRenderer;

public class VectorsRenderer extends ShapeRenderer {

  private final static float arrowHeadOffset = -0.2f;
  private final Vector3f vector2 = new Vector3f();
  private final Point3f pointVectorEnd = new Point3f();
  private final Point3f pointArrowHead = new Point3f();
  private final Point3i screenVectorEnd = new Point3i();
  private final Point3i screenArrowHead = new Point3i();
  private final Vector3f headOffsetVector = new Vector3f();
  
  private int diameter;
  //float headWidthAngstroms;
  private int headWidthPixels;
  private float vectorScale;
  private boolean vectorSymmetry;
  private float headScale;
  private boolean doShaft;


  @Override
  protected void render() {
    Vectors vectors = (Vectors) shape;
    if (!vectors.isActive)
      return;
    short[] mads = vectors.mads;
    if (mads == null)
      return;
    Atom[] atoms = vectors.atoms;
    short[] colixes = vectors.colixes;
    for (int i = modelSet.getAtomCount(); --i >= 0;) {
      Atom atom = atoms[i];
      if (!atom.isVisible(myVisibilityFlag))
        continue;
      Vector3f vibrationVector = viewer.getVibrationVector(i);
      if (vibrationVector == null)
        continue;
      vectorScale = viewer.getVectorScale();
      vectorSymmetry = viewer.getVectorSymmetry();
      if (transform(mads[i], atom, vibrationVector)
          && g3d.setColix(Shape.getColix(colixes, i, atom))) {
        renderVector(atom);
        if (vectorSymmetry) {
          vector2.set(vibrationVector);
          vector2.scale(-1);
          transform(mads[i], atom, vector2);
          renderVector(atom);          
        }
      }
    }
  }

  private boolean transform(short mad, Atom atom, Vector3f vibrationVector) {
    float len = vibrationVector.length();
    // to have the vectors move when vibration is turned on
    if (Math.abs(len * vectorScale) < 0.01)
      return false;
    headScale = arrowHeadOffset;
    if (vectorScale < 0)
      headScale = -headScale;
    doShaft = (0.1 + Math.abs(headScale/len) < Math.abs(vectorScale));
    headOffsetVector.set(vibrationVector);
    headOffsetVector.scale(headScale / len);
    pointVectorEnd.scaleAdd(vectorScale, vibrationVector, atom);
    pointArrowHead.set(pointVectorEnd);
    pointArrowHead.add(headOffsetVector);
    screenArrowHead.set(viewer.transformPoint(pointArrowHead, vibrationVector));
    screenVectorEnd.set(viewer.transformPoint(pointVectorEnd, vibrationVector));
    diameter = (mad < 1 ? 1 : mad <= 20 ? mad : viewer.scaleToScreen(screenVectorEnd.z, mad));
    headWidthPixels = (int)(diameter * 2.0f);
    if (headWidthPixels < diameter + 2)
      headWidthPixels = diameter + 2;
    return true;
  }
  
  private void renderVector(Atom atom) {
    if (doShaft)
      g3d.fillCylinderScreen(Graphics3D.ENDCAPS_OPEN, diameter, atom.screenX,
          atom.screenY, atom.screenZ, screenArrowHead.x, screenArrowHead.y,
          screenArrowHead.z);
    g3d.fillConeScreen(Graphics3D.ENDCAPS_FLAT, headWidthPixels, screenArrowHead,
        screenVectorEnd, false);
  }
}
