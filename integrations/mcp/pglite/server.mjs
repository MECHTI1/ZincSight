// One ephemeral PostgreSQL instance per prediction, accessible only on loopback.
import { readFile, writeFile, rename } from 'node:fs/promises';
import { PGlite } from '@electric-sql/pglite';
import { PGLiteSocketServer } from '@electric-sql/pglite-socket';

const [dumpPath, readyPath] = process.argv.slice(2);
if (!dumpPath || !readyPath) throw new Error('Expected bundled SQL dump and readiness file paths.');
const db = await PGlite.create();

// Load ZincSight's bundled pg_dump text. PGlite accepts COPY through /dev/blob;
// psql's interactive COPY FROM stdin protocol is not supported by this bridge.
// This parser intentionally handles only the repository's four-table fixture.
const sql = await readFile(dumpPath, 'utf8');
const copies = /^COPY (.+) FROM stdin;\r?\n([\s\S]*?)^\\\.\r?$/gm;
let end = 0;
let count = 0;
for (const match of sql.matchAll(copies)) {
  await db.exec(sql.slice(end, match.index));
  await db.query(`COPY ${match[1]} FROM '/dev/blob';`, [], {blob: new Blob([match[2]])});
  end = match.index + match[0].length;
  count++;
}
if (count !== 4) throw new Error(`Expected four bundled COPY blocks, found ${count}.`);
await db.exec(sql.slice(end));

const server = new PGLiteSocketServer({db, host: '127.0.0.1', port: 0, maxConnections: 64});
await server.start();
const version = (await db.query('SELECT version()')).rows[0].version;
const uri = `postgresql://postgres@${server.getServerConn()}/postgres?sslmode=disable`;
await writeFile(`${readyPath}.tmp`, JSON.stringify({uri, version}));
await rename(`${readyPath}.tmp`, readyPath);

for (const signal of ['SIGINT', 'SIGTERM']) process.on(signal, async () => {
  await server.stop();
  await db.close();
  process.exit(0);
});
