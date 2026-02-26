import { useRef, useMemo } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, Grid } from '@react-three/drei';
import * as THREE from 'three';

// ---------------------------------------------------------------------------
// Props
// ---------------------------------------------------------------------------

interface TurbineViewer3DProps {
  towerHeight?: number;
  rotorDiameter?: number;
  hubHeight?: number;
  nacelleLength?: number;
  className?: string;
}

// ---------------------------------------------------------------------------
// Tower geometry: tapered cylinder from base to top
// ---------------------------------------------------------------------------

function Tower({
  height,
  baseRadius,
  topRadius,
}: {
  height: number;
  baseRadius: number;
  topRadius: number;
}) {
  const geometry = useMemo(() => {
    return new THREE.CylinderGeometry(topRadius, baseRadius, height, 24, 1);
  }, [height, baseRadius, topRadius]);

  return (
    <mesh geometry={geometry} position={[0, height / 2, 0]} castShadow receiveShadow>
      <meshStandardMaterial color="#b0bec5" metalness={0.6} roughness={0.3} />
    </mesh>
  );
}

// ---------------------------------------------------------------------------
// Nacelle: box on top of the tower
// ---------------------------------------------------------------------------

function Nacelle({
  position,
  length,
}: {
  position: [number, number, number];
  length: number;
}) {
  const width = length * 0.35;
  const height = length * 0.35;

  return (
    <mesh position={[position[0] + length * 0.15, position[1], position[2]]} castShadow>
      <boxGeometry args={[length, height, width]} />
      <meshStandardMaterial color="#78909c" metalness={0.5} roughness={0.4} />
    </mesh>
  );
}

// ---------------------------------------------------------------------------
// Single Blade: flat rectangle with twist
// ---------------------------------------------------------------------------

function Blade({ length, angle }: { length: number; angle: number }) {
  const meshRef = useRef<THREE.Mesh>(null);

  const geometry = useMemo(() => {
    const bladeWidth = length * 0.08;
    const segments = 16;
    const geo = new THREE.BufferGeometry();

    const positions: number[] = [];
    const normals: number[] = [];
    const indices: number[] = [];

    for (let i = 0; i <= segments; i++) {
      const t = i / segments;
      const y = t * length;
      // Chord tapers linearly from root to tip
      const chord = bladeWidth * (1 - 0.7 * t);
      // Twist from 20 deg at root to 0 deg at tip
      const twist = ((1 - t) * 20 * Math.PI) / 180;

      const halfChord = chord / 2;

      // Left edge
      const lx = -halfChord * Math.cos(twist);
      const lz = -halfChord * Math.sin(twist);
      // Right edge
      const rx = halfChord * Math.cos(twist);
      const rz = halfChord * Math.sin(twist);

      positions.push(lx, y, lz);
      positions.push(rx, y, rz);

      // Normal pointing roughly outward
      normals.push(0, 0, 1);
      normals.push(0, 0, 1);

      if (i < segments) {
        const base = i * 2;
        indices.push(base, base + 1, base + 2);
        indices.push(base + 1, base + 3, base + 2);
      }
    }

    geo.setIndex(indices);
    geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geo.setAttribute('normal', new THREE.Float32BufferAttribute(normals, 3));
    geo.computeVertexNormals();

    return geo;
  }, [length]);

  return (
    <group rotation={[0, 0, (angle * Math.PI) / 180]}>
      <mesh ref={meshRef} geometry={geometry} castShadow>
        <meshStandardMaterial
          color="#eceff1"
          metalness={0.3}
          roughness={0.5}
          side={THREE.DoubleSide}
        />
      </mesh>
    </group>
  );
}

// ---------------------------------------------------------------------------
// Rotor: hub + 3 blades with rotation animation
// ---------------------------------------------------------------------------

function Rotor({
  position,
  bladeLength,
}: {
  position: [number, number, number];
  bladeLength: number;
}) {
  const groupRef = useRef<THREE.Group>(null);

  useFrame((_, delta) => {
    if (groupRef.current) {
      groupRef.current.rotation.x += delta * 0.8;
    }
  });

  const hubRadius = bladeLength * 0.04;

  return (
    <group position={position}>
      {/* Hub */}
      <mesh castShadow>
        <sphereGeometry args={[hubRadius, 16, 16]} />
        <meshStandardMaterial color="#90a4ae" metalness={0.6} roughness={0.3} />
      </mesh>

      {/* Blade rotor group - rotates around X axis (wind direction) */}
      <group ref={groupRef} rotation={[Math.PI / 2, 0, 0]}>
        {[0, 120, 240].map((angle) => (
          <Blade key={angle} length={bladeLength} angle={angle} />
        ))}
      </group>
    </group>
  );
}

// ---------------------------------------------------------------------------
// Scene
// ---------------------------------------------------------------------------

function TurbineScene({
  towerHeight,
  rotorDiameter,
  hubHeight,
  nacelleLength,
}: Required<Omit<TurbineViewer3DProps, 'className'>>) {
  const bladeLength = rotorDiameter / 2;
  const towerBaseRadius = towerHeight * 0.035;
  const towerTopRadius = towerHeight * 0.02;

  // Position the nacelle at hub height
  const nacelleY = hubHeight;
  const hubX = nacelleLength * 0.5;

  return (
    <>
      {/* Lighting */}
      <ambientLight intensity={0.4} />
      <directionalLight
        position={[50, 80, 30]}
        intensity={1}
        castShadow
        shadow-mapSize-width={1024}
        shadow-mapSize-height={1024}
      />
      <directionalLight position={[-30, 40, -20]} intensity={0.3} />

      {/* Tower */}
      <Tower
        height={towerHeight}
        baseRadius={towerBaseRadius}
        topRadius={towerTopRadius}
      />

      {/* Nacelle */}
      <Nacelle position={[0, nacelleY, 0]} length={nacelleLength} />

      {/* Rotor */}
      <Rotor
        position={[hubX + nacelleLength * 0.15, nacelleY, 0]}
        bladeLength={bladeLength}
      />

      {/* Ground grid */}
      <Grid
        args={[300, 300]}
        cellSize={10}
        cellThickness={0.5}
        cellColor="#334155"
        sectionSize={50}
        sectionThickness={1}
        sectionColor="#475569"
        fadeDistance={200}
        fadeStrength={1}
        position={[0, 0, 0]}
      />

      {/* Controls */}
      <OrbitControls
        enablePan
        enableZoom
        enableRotate
        target={[0, hubHeight * 0.5, 0]}
        minDistance={20}
        maxDistance={500}
        maxPolarAngle={Math.PI * 0.85}
      />
    </>
  );
}

// ---------------------------------------------------------------------------
// Exported component
// ---------------------------------------------------------------------------

export default function TurbineViewer3D({
  towerHeight = 87.6,
  rotorDiameter = 126,
  hubHeight = 90,
  nacelleLength = 12,
  className,
}: TurbineViewer3DProps) {
  return (
    <div className={className ?? 'h-full w-full min-h-[400px]'}>
      <Canvas
        camera={{
          position: [150, 80, 100],
          fov: 45,
          near: 1,
          far: 1000,
        }}
        shadows
        gl={{ antialias: true, alpha: false }}
        style={{ background: '#0f172a' }}
      >
        <TurbineScene
          towerHeight={towerHeight}
          rotorDiameter={rotorDiameter}
          hubHeight={hubHeight}
          nacelleLength={nacelleLength}
        />
      </Canvas>
    </div>
  );
}
